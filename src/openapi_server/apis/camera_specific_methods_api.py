# Coding: utf-8

from datetime import UTC, datetime

import astropy.units as u
import cabaret
import numpy as np
import yaml
from alpaca.exceptions import *
from alpaca.focuser import *
from alpaca.telescope import *
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time
from fastapi import APIRouter, Form, HTTPException, Path, Query  # Noqa: F401
from fastapi.responses import StreamingResponse

from openapi_server.models.bool_response import BoolResponse
from openapi_server.models.double_response import DoubleResponse
from openapi_server.models.image_array_response import ImageArrayResponse
from openapi_server.models.int_response import IntResponse
from openapi_server.models.method_response import MethodResponse
from openapi_server.models.string_array_response import StringArrayResponse
from openapi_server.models.string_response import StringResponse

# load config file
with open("config.yml", "r") as stream:
    config = yaml.safe_load(stream)

telescope = Telescope(config["TELESCOPE_IP"], config["TELESCOPE_DEVICE_NUMBER"])
focuser = Focuser(config["FOCUSER_IP"], config["FOCUSER_DEVICE_NUMBER"])

cabaret_camera = cabaret.Camera(
    name=config["CAMERA"]["name"],
    width=config["CAMERA"]["width"],  # pixels
    height=config["CAMERA"]["height"],  # pixels
    bin_x=config["CAMERA"]["bin_x"],  # binning factor in x
    bin_y=config["CAMERA"]["bin_y"],  # binning factor in y
    pitch=config["CAMERA"]["pitch"],  # pixel pitch, microns
    max_adu=config["CAMERA"]["max_adu"],  # maximum ADU value
    well_depth=config["CAMERA"]["well_depth"],  # electrons
    bias=config["CAMERA"]["bias"],  # ADU
    gain=config["CAMERA"]["gain"],  # e-/ADU
    read_noise=config["CAMERA"]["read_noise"],  # e-
    dark_current=config["CAMERA"]["dark_current"],  # e-/s
    average_quantum_efficiency=config["CAMERA"][
        "average_quantum_efficiency"
    ],  # fraction
)

cabaret_site = cabaret.Site(
    sky_background=config["SITE"]["sky_background"],  # e-/m2/arcsec2/s
    seeing=config["SITE"]["seeing"],  # arcsec
)

cabaret_telescope = cabaret.Telescope(
    focal_length=telescope.FocalLength,  # meters
    diameter=telescope.ApertureDiameter,  # meters
)

print(cabaret_telescope)


class ASCOM_Camera:
    device_number: int = 0
    width: int = cabaret_camera.width
    height: int = cabaret_camera.height
    startx: int = 0
    starty: int = 0
    numx: int = cabaret_camera.width
    numy: int = cabaret_camera.height
    binx: int = cabaret_camera.bin_x
    biny: int = cabaret_camera.bin_y
    gain: float = cabaret_camera.gain
    pitch: float = cabaret_camera.pitch
    welldepth: int = cabaret_camera.well_depth
    bias: int = cabaret_camera.bias
    maxadu: int = cabaret_camera.max_adu
    sensor_name: str = cabaret_camera.name
    imageready: bool = False
    imagearray: np.ndarray = None
    lastexposurestarttime: datetime = None
    lastexposureduration: float = None
    state: int = 0
    temperature: float = -10
    cooleron: bool = False


camera = ASCOM_Camera()


# @dataclass
# class Camera:
#     name: str = "gaia-camera-simulated"
#     width: int = 1024  # pixels
#     height: int = 1024  # pixels
#     bin_x: int = 1  # binning factor in x
#     bin_y: int = 1  # binning factor in y
#     pitch: float = 13.5  # pixel pitch, microns
#     plate_scale: float | None = (
#         None  # arcsec/pixel (calculated from pitch+telescope if None)
#     )
#     max_adu: int = 2**16 - 1  # maximum ADU value
#     well_depth: int = 2**16 - 1  # electrons
#     bias: int = 300  # ADU
#     gain: float = 1.0  # e-/ADU
#     read_noise: float = 6.2  # e-
#     dark_current: float = 0.2  # e-/s
#     average_quantum_efficiency: float = 0.8  # fraction


# @dataclass
# class Telescope:
#     focal_length: float = 8.0  # meters
#     diameter: float = 1.0  # meters
#     collecting_area: float | None = None  # m2 (calculated from diameter if None)


# @dataclass
# class Site:
#     sky_background: float = 150  # for I+z band in Paranal, e-/m2/arcsec2/s
#     seeing: float = 1.3  # arcsec


print("Starting simulator with config:")
print(config)


if config["FLATS"]:
    obs_lat = telescope.SiteLatitude
    obs_lon = telescope.SiteLongitude
    obs_alt = telescope.SiteElevation

start_time = datetime.now(UTC)


def tracking_error():
    ## add 1 arcsec drift in RA and 0.5 arcsec drift in Dec every second
    t = (datetime.now(UTC) - start_time).total_seconds()

    # sum up to get total error
    return 1 / 3600 * t, 0.5 / 3600 * t


def generate_image(exp_time, light=1):

    if light == 1:

        # get telescope position
        ra = (telescope.RightAscension / 24) * 360
        dec = telescope.Declination

        seeing_multiplier = (
            1
            + np.abs(focuser.Position - config["FOCUSER"]["sharp_pos"])
            / config["FOCUSER"]["blurred_offset"]
        )

        if seeing_multiplier > 5:
            seeing_multiplier = 5

        cabaret_site.seeing = config["SITE"]["seeing"] * seeing_multiplier  # arcsec

        print()
        print("Telescope position:")
        print(ra, dec)

        # add tracking error
        if config["TRACKING_ERROR"]:
            tracking_error_ra, tracking_error_dec = tracking_error()

            ra += tracking_error_ra
            dec += tracking_error_dec

        image = cabaret.generate_image(
            ra=ra,
            dec=dec,
            exp_time=exp_time,
            dateobs=datetime.now(UTC),
            light=light,
            camera=cabaret_camera,
            site=cabaret_site,
            telescope=cabaret_telescope,
            tmass=config["TMASS"],
        )

        if config["FLATS"]:

            # add sun scattered light during twilight for flats
            obs_time = Time.now()  # Time('2022-08-06T22:52:57.000', scale='utc')
            sun_position = get_sun(obs_time)
            obs_location = EarthLocation(
                lat=obs_lat * u.deg, lon=obs_lon * u.deg, height=obs_alt * u.m
            )
            sun_altaz = sun_position.transform_to(
                AltAz(obstime=obs_time, location=obs_location)
            )
            sun_alt = sun_altaz.alt.deg

            # from https://arxiv.org/pdf/1407.8283.pdf scaled to match a zYJ obs (~calibrated from single flat. This is fine for testing.)
            flux_sky = (
                (cabaret_camera.plate_scale) ** 2
                * 100
                * 10 ** (0.415 * sun_alt + 5.926)
                * cabaret_camera.average_quantum_efficiency
                * cabaret_telescope.collecting_area
                * exp_time
            )  # [electrons pixel^-1 ]

            if flux_sky > 1e9:
                flux_sky = 1e9

            # add sky background with poisson noise
            image += np.random.poisson(flux_sky, (camera.numy, camera.numx)).astype(
                np.uint16
            )

    else:
        # make base image with only dark current
        image = cabaret.generate_image(
            ra=0,
            dec=0,
            exp_time=exp_time,
            dateobs=datetime.now(UTC),
            light=light,
            camera=cabaret_camera,
            site=cabaret_site,
            telescope=cabaret_telescope,
        )

    # clip image to max adu
    image = np.clip(image, 0, cabaret_camera.max_adu)

    # save as fits
    # if config['SAVE_IMAGES']:
    #     hdu = fits.PrimaryHDU(image)
    #     hdr = hdu.header
    #     hdr.update(wcs.to_header())
    #     hdr['EXPTIME'] = exp_time
    #     hdr['DATE-OBS'] = datetime.now(UTC).isoformat()
    #     hdr['INSTRUME'] = camera.sensor_name
    #     hdr['TELESCOP'] = 'Gaia'
    #     hdr['RA'] = center.ra.deg
    #     hdr['DEC'] = center.dec.deg
    #     hdu.writeto('latest.fits', overwrite=True)

    return image


async def bytes_generator(image_array):

    b = int(1).to_bytes(4, "little")  # metaversion
    b += int(0).to_bytes(4, "little")  # error
    b += int(0).to_bytes(4, "little")  # clientid
    b += int(0).to_bytes(4, "little")  # deviceid
    b += int(44).to_bytes(4, "little")  # DataStart
    b += int(2).to_bytes(4, "little")  # ImageElementType
    b += int(8).to_bytes(4, "little")  # TransmissionElementType
    b += int(2).to_bytes(4, "little")  # Rank
    b += int(camera.numx).to_bytes(4, "little")  # Dimension1
    b += int(camera.numy).to_bytes(4, "little")  # Dimension2
    b += int(0).to_bytes(4, "little")  # Dimension3

    b += image_array.T.tobytes()

    yield b


router = APIRouter()
### USED METHODS ###


@router.put(
    "/api/v1/camera/{device_number}/abortexposure",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Aborts the current exposure",
    response_model_by_alias=True,
)
async def abortexposure_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """Aborts the current exposure, if any, and returns the camera to IDle state."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    camera.state = 0

    return MethodResponse(
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/binx",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the binning factor for the X axis.",
    response_model_by_alias=True,
)
async def binx_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the binning factor for the X axis."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return IntResponse(
        Value=camera.binx,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.put(
    "/api/v1/camera/{device_number}/binx",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Sets the binning factor for the X axis.",
    response_model_by_alias=True,
)
async def binx_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    BinX: int = Form(1, description="The X binning value"),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """Sets the binning factor for the X axis."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    camera.binx = BinX

    return MethodResponse(
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/biny",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the binning factor for the Y axis.",
    response_model_by_alias=True,
)
async def biny_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the binning factor for the Y axis."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return IntResponse(
        Value=camera.biny,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.put(
    "/api/v1/camera/{device_number}/biny",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Sets the binning factor for the Y axis.",
    response_model_by_alias=True,
)
async def biny_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    BinY: int = Form(1, description="The Y binning value"),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """Sets the binning factor for the Y axis."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    camera.biny = BinY

    return MethodResponse(
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/camerastate",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the camera operational state.",
    response_model_by_alias=True,
)
async def camerastate_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the current camera operational state as an integer. 0 &#x3D; CameraIDle , 1 &#x3D; CameraWaiting , 2 &#x3D; CameraExposing , 3 &#x3D; CameraReading , 4 &#x3D; CameraDownload , 5 &#x3D; CameraError"""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return IntResponse(
        Value=camera.state,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/cameraxsize",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the width of the CCD camera chip.",
    response_model_by_alias=True,
)
async def cameraxsize_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the width of the CCD camera chip in unbinned pixels."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return IntResponse(
        Value=camera.width,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/cameraysize",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the height of the CCD camera chip.",
    response_model_by_alias=True,
)
async def cameraysize_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the height of the CCD camera chip in unbinned pixels."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return IntResponse(
        Value=camera.height,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/ccdtemperature",
    responses={
        200: {"model": DoubleResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the current CCD temperature",
    response_model_by_alias=True,
)
async def ccdtemperature_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> DoubleResponse:
    """Returns the current CCD temperature in degrees Celsius."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return DoubleResponse(
        Value=camera.temperature,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/cooleron",
    responses={
        200: {"model": BoolResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the current cooler on/off state.",
    response_model_by_alias=True,
)
async def cooleron_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> BoolResponse:
    """Returns the current cooler on/off state."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return BoolResponse(
        Value=camera.cooleron,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.put(
    "/api/v1/camera/{device_number}/cooleron",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Turns the camera cooler on and off",
    response_model_by_alias=True,
)
async def cooleron_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    CoolerOn: bool = Form(True, description="Cooler state"),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """Turns on and off the camera cooler. True &#x3D; cooler on, False &#x3D; cooler off"""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    camera.cooleron = CoolerOn

    return MethodResponse(
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/imagearray",
    responses={
        200: {"model": bytes, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns an array of integers containing the exposure pixel values",
    response_model_by_alias=True,
)
async def imagearray_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> StreamingResponse:
    """Returns an array of 32bit integers containing the pixel values from the last exposure. This call can return either a 2 dimension (monochrome images) or 3 dimension (colour or multi-plane images) array of size NumX \\* NumY or NumX \\* NumY \\* NumPlanes. Where applicable, the size of NumPlanes has to be determined by inspection of the returned Array.  Since 32bit integers are always returned by this call, the returned JSON Type value (0 &#x3D; Unknown, 1 &#x3D; short(16bit), 2 &#x3D; int(32bit), 3 &#x3D; Double) is always 2. The number of planes is given in the returned Rank value.  When de-serialising to an object it is essential to know the array Rank beforehand so that the correct data class can be used. This can be achieved through a regular expression or by direct parsing of the returned JSON string to extract the Type and Rank values before de-serialising.  This regular expression accomplishes the extraction into two named groups Type and Rank, which can then be used to select the correct de-serialisation data Class:  __&#x60;^*\&quot;Type\&quot;:(?&lt;Type&gt;\\d*),\&quot;Rank\&quot;:(?&lt;Rank&gt;\\d*)&#x60;__  When the SensorType is Monochrome, RGGB, CMYG, CMYG2 or LRGB, the serialised JSON array should have 2 dimensions. For example, the returned array should appear as below if NumX &#x3D; 7, NumY &#x3D; 5  and Pxy represents the pixel value at the zero based position x across and y down the image with the origin in the top left corner of the image.    Please note that this is \&quot;column-major\&quot; order (column changes most rapidly) from the image's row and column perspective, while, from the array's perspective, serialisation is actually effected in \&quot;row-major\&quot; order (rightmost index changes most rapidly).  This unintuitive outcome arises because the ASCOM Camera Interface specification defines the image column dimension as the rightmost array dimension.  [  [P00,P01,P02,P03,P04],  [P10,P11,P12,P13,P14],  [P20,P21,P22,P23,P24],  [P30,P31,P32,P33,P34],  [P40,P41,P42,P43,P44],  [P50,P51,P52,P53,P54],  [P60,P61,P62,P63,P64]  ]  When the SensorType is Color, the serialised JSON array will have 3 dimensions. For example, the returned array should appear as below if NumX &#x3D; 7, NumY &#x3D; 5  and Rxy, Gxy and Bxy represent the red, green and blue pixel values at the zero based position x across and y down the image with the origin in the top left corner of the image.  Please see note above regarding element ordering.  [  [[R00,G00,B00],[R01,G01,B01],[R02,G02,B02],[R03,G03,B03],[R04,G04,B04]],  [[R10,G10,B10],[R11,G11,B11],[R12,G12,B12],[R13,G13,B13],[R14,G14,B14]],  [[R20,G20,B20],[R21,G21,B21],[R22,G22,B22],[R23,G23,B23],[R24,G24,B24]],  [[R30,G30,B30],[R31,G31,B31],[R32,G32,B32],[R33,G33,B33],[R34,G34,B34]],  [[R40,G40,B40],[R41,G41,B41],[R42,G42,B42],[R43,G43,B43],[R44,G44,B44]],  [[R50,G50,B50],[R51,G51,B51],[R52,G52,B52],[R53,G53,B53],[R54,G54,B54]],  [[R60,G60,B60],[R61,G61,B61],[R62,G62,B62],[R63,G63,B63],[R64,G64,B64]],  ]  __&#x60;Performance&#x60;__  Returning an image from an Alpaca device as a JSON array is very inefficient and can result in delays of 30 or more seconds while client and device process and send the huge JSON string over the network.  A new, much faster mechanic called ImageBytes - [Alpaca ImageBytes Concepts and Implementation](https://www.ascom-standards.org/Developer/AlpacaImageBytes.pdf) has been developed that sends data as a binary byte stream and can offer a 10 to 20 fold reduction in transfer time.  It is strongly recommended that Alpaca Cameras implement the ImageBytes mechanic as well as the JSON mechanic."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    headers = {"Content-Type": "application/imagebytes"}
    return StreamingResponse(
        bytes_generator(camera.imagearray),
        headers=headers,
        media_type="application/imagebytes",
    )


@router.get(
    "/api/v1/camera/{device_number}/imageready",
    responses={
        200: {"model": BoolResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Indicates that an image is ready to be downloaded",
    response_model_by_alias=True,
)
async def imageready_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> BoolResponse:
    """Returns a flag indicating whether the image is ready to be downloaded from the camera."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    if (
        datetime.now(UTC) - camera.lastexposurestarttime
    ).total_seconds() > camera.lastexposureduration:
        camera.imageready = True
        camera.state = 0
    else:
        camera.imageready = False

    return BoolResponse(
        Value=camera.imageready,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/lastexposurestarttime",
    responses={
        200: {"model": StringResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Start time of the last exposure in FITS standard format.",
    response_model_by_alias=True,
)
async def lastexposurestarttime_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> StringResponse:
    """Reports the actual exposure start in the FITS-standard CCYY-MM-DDThh:mm:ss[.sssif camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")] format.
    """
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return StringResponse(
        Value=str(camera.lastexposurestarttime),
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/numx",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the current subframe width",
    response_model_by_alias=True,
)
async def numx_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the current subframe width, if binning is active, value is in binned pixels."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return IntResponse(
        Value=camera.numx,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.put(
    "/api/v1/camera/{device_number}/numx",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Sets the current subframe width",
    response_model_by_alias=True,
)
async def numx_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    NumX: int = Form(
        0,
        description="Sets the subframe width, if binning is active, value is in binned pixels.",
    ),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """Sets the current subframe width."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    camera.numx = NumX

    return MethodResponse(
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/numy",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the current subframe height",
    response_model_by_alias=True,
)
async def numy_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the current subframe height, if binning is active, value is in binned pixels."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return IntResponse(
        Value=camera.numy,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.put(
    "/api/v1/camera/{device_number}/numy",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Sets the current subframe height",
    response_model_by_alias=True,
)
async def numy_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    NumY: int = Form(
        0,
        description="Sets the subframe height, if binning is active, value is in binned pixels.",
    ),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """Sets the current subframe height."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    camera.numy = NumY

    return MethodResponse(
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/pixelsizex",
    responses={
        200: {"model": DoubleResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Width of CCD chip pixels (microns)",
    response_model_by_alias=True,
)
async def pixelsizex_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> DoubleResponse:
    """Returns the width of the CCD chip pixels in microns."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return DoubleResponse(
        Value=camera.pitch,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/pixelsizey",
    responses={
        200: {"model": DoubleResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Height of CCD chip pixels (microns)",
    response_model_by_alias=True,
)
async def pixelsizey_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> DoubleResponse:
    """Returns the Height of the CCD chip pixels in microns."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return DoubleResponse(
        Value=camera.pitch,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/sensorname",
    responses={
        200: {"model": StringResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Sensor name",
    response_model_by_alias=True,
)
async def sensorname_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> StringResponse:
    """The name of the sensor used within the camera."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return StringResponse(
        Value=camera.sensor_name,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/sensortype",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Type of information returned by the the camera sensor (monochrome or colour)",
    response_model_by_alias=True,
)
async def sensortype_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns a value indicating whether the sensor is monochrome, or what Bayer matrix it encodes. Where:  - 0 &#x3D; Monochrome, - 1 &#x3D; Colour not requiring Bayer decoding - 2 &#x3D; RGGB Bayer encoding - 3 &#x3D; CMYG Bayer encoding - 4 &#x3D; CMYG2 Bayer encoding - 5 &#x3D; LRGB TRUESENSE Bayer encoding.  Please see the ASCOM Help fie for more informaiton on the SensorType."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return IntResponse(
        Value=0,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/setccdtemperature",
    responses={
        200: {"model": DoubleResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the current camera cooler setpoint in degrees Celsius.",
    response_model_by_alias=True,
)
async def setccdtemperature_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> DoubleResponse:
    """Returns the current camera cooler setpoint in degrees Celsius."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return DoubleResponse(
        Value=camera.temperature,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/cansetccdtemperature",
    responses={
        200: {"model": BoolResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns a flag indicating whether this camera supports setting the CCD temperature",
    response_model_by_alias=True,
)
async def cansetccdtemperature_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> BoolResponse:
    """Returns a flag indicatig whether this camera supports setting the CCD temperature"""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return BoolResponse(
        Value=True,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.put(
    "/api/v1/camera/{device_number}/setccdtemperature",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Set the camera's cooler setpoint (degrees Celsius).",
    response_model_by_alias=True,
)
async def setccdtemperature_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    SetCCDTemperature: float = Form(
        -10, description="Temperature set point (degrees Celsius)."
    ),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """Set's the camera's cooler setpoint in degrees Celsius."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    camera.temperature = SetCCDTemperature

    return MethodResponse(
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.put(
    "/api/v1/camera/{device_number}/startexposure",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Starts an exposure",
    response_model_by_alias=True,
)
async def startexposure_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    Duration: float = Form(None, description="Duration of exposure in seconds"),
    Light: bool = Form(None, description="True if light frame, false if dark frame."),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """Starts an exposure. Use ImageReady to check when the exposure is complete."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    camera.imageready = False

    camera.lastexposurestarttime = datetime.now(UTC)

    camera.lastexposureduration = Duration

    camera.state = 2

    camera.imagearray = generate_image(Duration, Light)

    return MethodResponse(
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/startx",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Return the current subframe X axis start position",
    response_model_by_alias=True,
)
async def startx_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Sets the subframe start position for the X axis (0 based) and returns the current value. If binning is active, value is in binned pixels."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return IntResponse(
        Value=camera.startx,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.put(
    "/api/v1/camera/{device_number}/startx",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Sets the current subframe X axis start position",
    response_model_by_alias=True,
)
async def startx_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    StartX: int = Form(
        0, description="The subframe X axis start position in binned pixels."
    ),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """Sets the current subframe X axis start position in binned pixels."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    camera.startx = StartX

    return MethodResponse(
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/starty",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Return the current subframe Y axis start position",
    response_model_by_alias=True,
)
async def starty_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Sets the subframe start position for the Y axis (0 based) and returns the current value. If binning is active, value is in binned pixels."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return IntResponse(
        Value=camera.starty,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.put(
    "/api/v1/camera/{device_number}/starty",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Sets the current subframe Y axis start position",
    response_model_by_alias=True,
)
async def starty_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    StartY: int = Form(
        0, description="The subframe Y axis start position in binned pixels."
    ),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """Sets the current subframe Y axis start position in binned pixels."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    camera.starty = StartY
    return MethodResponse(
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


### UNUSED METHODS


@router.get(
    "/api/v1/camera/{device_number}/bayeroffsetx",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the X offset of the Bayer matrix.",
    response_model_by_alias=True,
)
async def bayeroffsetx_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the X offset of the Bayer matrix, as defined in SensorType."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/bayeroffsety",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the Y offset of the Bayer matrix.",
    response_model_by_alias=True,
)
async def bayeroffsety_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the Y offset of the Bayer matrix, as defined in SensorType."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/canasymmetricbin",
    responses={
        200: {"model": BoolResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Indicates whether the camera supports asymmetric binning",
    response_model_by_alias=True,
)
async def canasymmetricbin_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> BoolResponse:
    """Returns a flag showing whether this camera supports asymmetric binning"""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return BoolResponse(
        Value=False,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/canfastreadout",
    responses={
        200: {"model": BoolResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Indicates whether the camera has a fast readout mode.",
    response_model_by_alias=True,
)
async def canfastreadout_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> BoolResponse:
    """Indicates whether the camera has a fast readout mode."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    return BoolResponse(
        Value=False,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/canabortexposure",
    responses={
        200: {"model": BoolResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Indicates whether the camera can abort exposures.",
    response_model_by_alias=True,
)
async def canabortexposure_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> BoolResponse:
    """Returns true if the camera can abort exposures; false if not."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return BoolResponse(
        Value=False,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/cangetcoolerpower",
    responses={
        200: {"model": BoolResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Indicates whether the camera's cooler power setting can be read.",
    response_model_by_alias=True,
)
async def cangetcoolerpower_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> BoolResponse:
    """If true, the camera's cooler power setting can be read."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return BoolResponse(
        Value=False,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/canpulseguide",
    responses={
        200: {"model": BoolResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns a flag indicating whether this camera supports pulse guiding",
    response_model_by_alias=True,
)
async def canpulseguide_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> BoolResponse:
    """Returns a flag indicating whether this camera supports pulse guiding."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return BoolResponse(
        Value=False,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/canstopexposure",
    responses={
        200: {"model": BoolResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns a flag indicating whether this camera can stop an exposure that is in progress",
    response_model_by_alias=True,
)
async def canstopexposure_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> BoolResponse:
    """Returns a flag indicating whether this camera can stop an exposure that is in progress"""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return BoolResponse(
        Value=False,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/exposuremax",
    responses={
        200: {"model": DoubleResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the maximum exposure time supported by StartExposure.",
    response_model_by_alias=True,
)
async def exposuremax_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> DoubleResponse:
    """Returns the maximum exposure time supported by StartExposure."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/exposuremin",
    responses={
        200: {"model": DoubleResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the Minimium exposure time",
    response_model_by_alias=True,
)
async def exposuremin_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> DoubleResponse:
    """Returns the Minimium exposure time in seconds that the camera supports through StartExposure."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/electronsperadu",
    responses={
        200: {"model": DoubleResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the gain of the camera",
    response_model_by_alias=True,
)
async def electronsperadu_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> DoubleResponse:
    """Returns the gain of the camera in photoelectrons per A/D unit."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return DoubleResponse(
        Value=camera.gain,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/exposureresolution",
    responses={
        200: {"model": DoubleResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the smallest increment in exposure time supported by StartExposure.",
    response_model_by_alias=True,
)
async def exposureresolution_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> DoubleResponse:
    """Returns the smallest increment in exposure time supported by StartExposure."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/fastreadout",
    responses={
        200: {"model": BoolResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns whenther Fast Readout Mode is enabled.",
    response_model_by_alias=True,
)
async def fastreadout_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> BoolResponse:
    """Returns whenther Fast Readout Mode is enabled."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.put(
    "/api/v1/camera/{device_number}/fastreadout",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Sets whether Fast Readout Mode is enabled.",
    response_model_by_alias=True,
)
async def fastreadout_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    FastReadout: bool = Form(False, description="True to enable fast readout mode"),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """Sets whether Fast Readout Mode is enabled."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/fullwellcapacity",
    responses={
        200: {"model": DoubleResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Reports the full well capacity of the camera",
    response_model_by_alias=True,
)
async def fullwellcapacity_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> DoubleResponse:
    """Reports the full well capacity of the camera in electrons, at the current camera settings (binning, SetupDialog settings, etc.)."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return DoubleResponse(
        Value=camera.welldepth,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/gain",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the camera's gain",
    response_model_by_alias=True,
)
async def gain_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """The camera's gain (GAIN VALUE MODE) OR the index of the selected camera gain description in the Gains array (GAINS INDEX MODE)."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.put(
    "/api/v1/camera/{device_number}/gain",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Sets the camera's gain.",
    response_model_by_alias=True,
)
async def gain_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    Gain: int = Form(
        0, description="Index of the current camera gain in the Gains string array."
    ),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """The camera's gain (GAIN VALUE MODE) OR the index of the selected camera gain description in the Gains array (GAINS INDEX MODE)."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/gainmax",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Maximum Gain value of that this camera supports",
    response_model_by_alias=True,
)
async def gainmax_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the maximum value of Gain."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/gainmin",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Minimum Gain value of that this camera supports",
    response_model_by_alias=True,
)
async def gainmin_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the Minimum value of Gain."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/gains",
    responses={
        200: {"model": StringArrayResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="List of Gain names supported by the camera",
    response_model_by_alias=True,
)
async def gains_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> StringArrayResponse:
    """Returns the Gains supported by the camera."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/hasshutter",
    responses={
        200: {"model": BoolResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Indicates whether the camera has a mechanical shutter",
    response_model_by_alias=True,
)
async def hasshutter_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> BoolResponse:
    """Returns a flag indicating whether this camera has a mechanical shutter."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return BoolResponse(
        Value=True,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/heatsinktemperature",
    responses={
        200: {"model": DoubleResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the current heat sink temperature.",
    response_model_by_alias=True,
)
async def heatsinktemperature_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> DoubleResponse:
    """Returns the current heat sink temperature (called \&quot;ambient temperature\&quot; by some manufacturers) in degrees Celsius."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/imagearrayvariant",
    responses={
        200: {"model": ImageArrayResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns an array of int containing the exposure pixel values",
    response_model_by_alias=True,
)
async def imagearrayvariant_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> ImageArrayResponse:
    """Returns an array containing the pixel values from the last exposure. This call can return either a 2 dimension (monochrome images) or 3 dimension (colour or multi-plane images) array of size NumX \\* NumY  or NumX \\* NumY \\* NumPlanes. Where applicable, the size of NumPlanes has to be determined by inspection of the returned Array.  This call can return values as short(16bit) integers, int(32bit) integers or double floating point values. The nature of the returned values is given in the Type Parameter: 0 &#x3D; Unknown, 1 &#x3D; short(16bit), 2 &#x3D; int(32bit), 3 &#x3D; Double. The number of planes is given in the returned Rank value.  When deserialising to an object it helps enormously to know the Type and Rank beforehand so that the correct data class can be used. This can be achieved through a regular expression or by direct parsing of the returned JSON string to extract the Type and Rank values before deserialising.  This regular expression accomplishes the extraction into two named groups Type and Rank, which can then be used to select the correct de-serialisation data Class:  __&#x60;^*\&quot;Type\&quot;:(?&lt;Type&gt;\\d*),\&quot;Rank\&quot;:(?&lt;Rank&gt;\\d*)&#x60;__  When the SensorType is Monochrome, RGGB, CMYG, CMYG2 or LRGB, the serialised JSON array should have 2 dimensions. For example, the returned array should appear as below if NumX &#x3D; 7, NumY &#x3D; 5  and Pxy represents the pixel value at the zero based position x across and y down the image with the origin in the top left corner of the image.    Please note that this is \&quot;column-major\&quot; order (column changes most rapidly) from the image's row and column perspective, while, from the array's perspective, serialisation is actually effected in \&quot;row-major\&quot; order (rightmost index changes most rapidly).  This unintuitive outcome arises because the ASCOM Camera Interface specification defines the image column dimension as the rightmost array dimension.  [  [P00,P01,P02,P03,P04],  [P10,P11,P12,P13,P14],  [P20,P21,P22,P23,P24],  [P30,P31,P32,P33,P34],  [P40,P41,P42,P43,P44],  [P50,P51,P52,P53,P54],  [P60,P61,P62,P63,P64]  ]  When the SensorType is Color, the serialised JSON array should have 3 dimensions. For example, the returned array should appear as below if NumX &#x3D; 7, NumY &#x3D; 5  and Rxy, Gxy and Bxy represent the red, green and blue pixel values at the zero based position x across and y down the image with the origin in the top left corner of the image.  Please see note above regarding element ordering.  [  [[R00,G00,B00],[R01,G01,B01],[R02,G02,B02],[R03,G03,B03],[R04,G04,B04]],  [[R10,G10,B10],[R11,G11,B11],[R12,G12,B12],[R13,G13,B13],[R14,G14,B14]],  [[R20,G20,B20],[R21,G21,B21],[R22,G22,B22],[R23,G23,B23],[R24,G24,B24]],  [[R30,G30,B30],[R31,G31,B31],[R32,G32,B32],[R33,G33,B33],[R34,G34,B34]],  [[R40,G40,B40],[R41,G41,B41],[R42,G42,B42],[R43,G43,B43],[R44,G44,B44]],  [[R50,G50,B50],[R51,G51,B51],[R52,G52,B52],[R53,G53,B53],[R54,G54,B54]],  [[R60,G60,B60],[R61,G61,B61],[R62,G62,B62],[R63,G63,B63],[R64,G64,B64]],  ]  __&#x60;Performance&#x60;__  Returning an image from an Alpaca device as a JSON array is very inefficient and can result in delays of 30 or more seconds while client and device process and send the huge JSON string over the network.  A new, much faster mechanic called ImageBytes - [Alpaca ImageBytes Concepts and Implementation](https://www.ascom-standards.org/Developer/AlpacaImageBytes.pdf) has been developed that sends data as a binary byte stream and can offer a 10 to 20 fold reduction in transfer time.  It is strongly recommended that Alpaca Cameras implement the ImageBytes mechanic as well as the JSON mechanic."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/ispulseguiding",
    responses={
        200: {"model": BoolResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Indicates that the camera is pulse guideing.",
    response_model_by_alias=True,
)
async def ispulseguiding_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> BoolResponse:
    """Returns a flag indicating whether the camera is currrently in a PulseGuide operation."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/maxadu",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Camera's maximum ADU value",
    response_model_by_alias=True,
)
async def maxadu_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Reports the maximum ADU value the camera can produce."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return IntResponse(
        Value=camera.maxadu,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/maxbinx",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Maximum  binning for the camera X axis",
    response_model_by_alias=True,
)
async def maxbinx_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the maximum allowed binning for the X camera axis"""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return IntResponse(
        Value=1,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/maxbiny",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Maximum  binning for the camera Y axis",
    response_model_by_alias=True,
)
async def maxbiny_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the maximum allowed binning for the Y camera axis"""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return IntResponse(
        Value=1,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/lastexposureduration",
    responses={
        200: {"model": DoubleResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Duration of the last exposure",
    response_model_by_alias=True,
)
async def lastexposureduration_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> DoubleResponse:
    """Reports the actual exposure duration in seconds (i.e. shutter open time)."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return DoubleResponse(
        Value=camera.lastexposureduration,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/offset",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the camera's offset",
    response_model_by_alias=True,
)
async def offset_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the camera's offset (OFFSET VALUE MODE) OR the index of the selected camera offset description in the offsets array (OFFSETS INDEX MODE)."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.put(
    "/api/v1/camera/{device_number}/offset",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Sets the camera's offset.",
    response_model_by_alias=True,
)
async def offset_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    Offset: int = Form(
        0, description="Index of the current camera offset in the offsets string array."
    ),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """Sets the camera's offset (OFFSET VALUE MODE) OR the index of the selected camera offset description in the offsets array (OFFSETS INDEX MODE)."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/offsetmax",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Maximum offset value of that this camera supports",
    response_model_by_alias=True,
)
async def offsetmax_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the maximum value of offset."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/offsetmin",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Minimum offset value of that this camera supports",
    response_model_by_alias=True,
)
async def offsetmin_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the Minimum value of offset."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/offsets",
    responses={
        200: {"model": StringArrayResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="List of offset names supported by the camera",
    response_model_by_alias=True,
)
async def offsets_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> StringArrayResponse:
    """Returns the offsets supported by the camera."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/percentcompleted",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Indicates percentage completeness of the current operation",
    response_model_by_alias=True,
)
async def percentcompleted_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """Returns the percentage of the current operation that is complete. If valid, returns an integer between 0 and 100, where 0 indicates 0% progress (function just started) and 100 indicates 100% progress (i.e. completion)."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.put(
    "/api/v1/camera/{device_number}/pulseguide",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Pulse guide in the specified direction for the specified time.",
    response_model_by_alias=True,
)
async def pulseguide_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    Direction: int = Form(
        None,
        description="Direction of movement (0 &#x3D; North, 1 &#x3D; South, 2 &#x3D; East, 3 &#x3D; West)",
    ),
    Duration: int = Form(None, description="Duration of movement in milli-seconds"),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """Activates the Camera's mount control sytem to instruct the mount to move in a particular direction for a given period of time"""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/readoutmode",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Indicates the canera's readout mode as an index into the array ReadoutModes",
    response_model_by_alias=True,
)
async def readoutmode_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> IntResponse:
    """ReadoutMode is an index into the array ReadoutModes and returns the desired readout mode for the camera. Defaults to 0 if not set."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return IntResponse(
        Value=0,
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.put(
    "/api/v1/camera/{device_number}/readoutmode",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Set the camera's readout mode",
    response_model_by_alias=True,
)
async def readoutmode_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ReadoutMode: int = Form(
        0,
        description="Index into the ReadoutModes array of string readout mode names indicating the camera's current readout mode.",
    ),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """Sets the ReadoutMode as an index into the array ReadoutModes."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")

    return MethodResponse(
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/readoutmodes",
    responses={
        200: {"model": StringArrayResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="List of available readout modes",
    response_model_by_alias=True,
)
async def readoutmodes_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> StringArrayResponse:
    """This property provides an array of strings, each of which describes an available readout mode of the camera. At least one string must be present in the list."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    return StringArrayResponse(
        Value=["no"],
        ClientTransactionID=ClientTransactionID,
        ServerTransactionID=ClientTransactionID,
        ErrorNumber=0,
        ErrorMessage="",
    )


@router.get(
    "/api/v1/camera/{device_number}/subexposureduration",
    responses={
        200: {"model": DoubleResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Camera's sub-exposure interval",
    response_model_by_alias=True,
)
async def subexposureduration_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> DoubleResponse:
    """The Camera's sub exposure duration in seconds. Only available in Camera Interface Version 3 and later."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.put(
    "/api/v1/camera/{device_number}/subexposureduration",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Sets the current Sub Exposure Duration",
    response_model_by_alias=True,
)
async def subexposureduration_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    SubExposureDuration: float = Form(
        0, description="The request sub exposure duration in seconds"
    ),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """Sets image sub exposure duration in seconds. Only available in Camera Interface Version 3 and later."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/camera/{device_number}/coolerpower",
    responses={
        200: {"model": DoubleResponse, "description": "Driver response"},
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Returns the present cooler power level",
    response_model_by_alias=True,
)
async def coolerpower_get(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Query(
        1,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
        ge=0,
        le=4294967295,
    ),
    ClientTransactionID: int = Query(
        1234,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
        ge=0,
        le=4294967295,
    ),
) -> DoubleResponse:
    """Returns the present cooler power level, in percent."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.put(
    "/api/v1/camera/{device_number}/stopexposure",
    responses={
        200: {
            "model": MethodResponse,
            "description": "Transaction complete or exception",
        },
        400: {
            "model": str,
            "description": "Method or parameter value error, check error message",
        },
        500: {
            "model": str,
            "description": "Server internal error, check error message",
        },
    },
    tags=["Camera Specific Methods"],
    summary="Stops the current exposure",
    response_model_by_alias=True,
)
async def stopexposure_put(
    device_number: int = Path(
        ...,
        description="Zero based device number as set on the server (0 to 4294967295)",
        ge=0,
        le=4294967295,
    ),
    ClientID: int = Form(
        None,
        description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.",
    ),
    ClientTransactionID: int = Form(
        None,
        description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.",
    ),
) -> MethodResponse:
    """Stops the current exposure, if any. If an exposure is in progress, the readout process is initiated. Ignored if readout is already in process."""
    if camera.device_number != device_number:
        raise HTTPException(status_code=404, detail="Device not found.")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")
