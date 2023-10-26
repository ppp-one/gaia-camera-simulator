# Coding: utf-8

from fastapi import (  # Noqa: F401
    APIRouter,
    Body,
    Form,
    Path,
    Query,
    HTTPException
)

from openapi_server.models.bool_response import BoolResponse
from openapi_server.models.device_type_device_number_commandblind_put_request import DeviceTypeDeviceNumberCommandblindPutRequest
from openapi_server.models.int_response import IntResponse
from openapi_server.models.method_response import MethodResponse
from openapi_server.models.string_array_response import StringArrayResponse
from openapi_server.models.string_response import StringResponse


router = APIRouter()


@router.put(
    "/api/v1/{device_type}/{device_number}/action",
    responses={
        200: {"model": StringResponse, "description": "Transaction complete or exception"},
        400: {"model": str, "description": "Method or parameter value error, check error message"},
        500: {"model": str, "description": "Server internal error, check error message"},
    },
    tags=["ASCOM Methods Common To All Devices"],
    summary="Invokes the named device-specific action.",
    response_model_by_alias=True,
)
async def device_type_device_number_action_put(
    device_type: str = Path(..., description="One of the recognised ASCOM device types e.g. telescope (must be lower case)", regex=r"^[a-z]*$"),
    device_number: int = Path(..., description="Zero based device number as set on the server (0 to 4294967295)", ge=0, le=4294967295),
    Action: str = Form(None, description="A well known name that represents the action to be carried out."),
    Parameters: str = Form(None, description="List of required parameters or an Empty String if none are required"),
    ClientID: int = Form(None, description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client."),
    ClientTransactionID: int = Form(None, description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs."),
) -> StringResponse:
    """Actions and SupportedActions are a standardised means for drivers to extend functionality beyond the built-in capabilities of the ASCOM device interfaces.  The key advantage of using Actions is that drivers can expose any device specific functionality required. The downside is that, in order to use these unique features, every application author would need to create bespoke code to present or exploit them.  The Action parameter and return strings are deceptively simple, but can support transmission of arbitrarily complex data structures, for example through JSON encoding.  This capability will be of primary value to  * &lt;span style&#x3D;\&quot;font-size:14px;\&quot;&gt;bespoke software and hardware configurations where a single entity controls both the consuming application software and the hardware / driver environment&lt;/span&gt;  * &lt;span style&#x3D;\&quot;font-size:14px;\&quot;&gt;a group of application and device authors to quickly formulate and try out new interface capabilities without requiring an immediate change to the ASCOM device interface, which will take a lot longer than just agreeing a name, input parameters and a standard response for an Action command.&lt;/span&gt;   The list of Action commands supported by a driver can be discovered through the SupportedActions property.  This method should return an error message and NotImplementedException error number (0x400) if the driver just implements the standard ASCOM device methods and has no bespoke, unique, functionality."""
    raise HTTPException(status_code=400, detail="Method not implemented.")


@router.put(
    "/api/v1/{device_type}/{device_number}/commandblind",
    responses={
        200: {"model": MethodResponse, "description": "Transaction complete or exception"},
        400: {"model": str, "description": "Method or parameter value error, check error message"},
        500: {"model": str, "description": "Server internal error, check error message"},
    },
    tags=["ASCOM Methods Common To All Devices"],
    summary="Transmits an arbitrary string to the device",
    response_model_by_alias=True,
)
async def device_type_device_number_commandblind_put(
    device_type: str = Path(..., description="One of the recognised ASCOM device types e.g. telescope (must be lower case)", regex=r"^[a-z]*$"),
    device_number: int = Path(..., description="Zero based device number as set on the server (0 to 4294967295)", ge=0, le=4294967295),
    device_type_device_number_commandblind_put_request: DeviceTypeDeviceNumberCommandblindPutRequest = Body(None, description=""),
) -> MethodResponse:
    """Transmits an arbitrary string to the device and does not wait for a response. Optionally, protocol framing characters may be added to the string before transmission."""
    raise HTTPException(status_code=400, detail="Method not implemented.")


@router.put(
    "/api/v1/{device_type}/{device_number}/commandbool",
    responses={
        200: {"model": BoolResponse, "description": "Transaction complete or exception"},
        400: {"model": str, "description": "Method or parameter value error, check error message"},
        500: {"model": str, "description": "Server internal error, check error message"},
    },
    tags=["ASCOM Methods Common To All Devices"],
    summary="Transmits an arbitrary string to the device and returns a boolean value from the device.",
    response_model_by_alias=True,
)
async def device_type_device_number_commandbool_put(
    device_type: str = Path(..., description="One of the recognised ASCOM device types e.g. telescope (must be lower case)", regex=r"^[a-z]*$"),
    device_number: int = Path(..., description="Zero based device number as set on the server (0 to 4294967295)", ge=0, le=4294967295),
    device_type_device_number_commandblind_put_request: DeviceTypeDeviceNumberCommandblindPutRequest = Body(None, description=""),
) -> BoolResponse:
    """Transmits an arbitrary string to the device and waits for a boolean response. Optionally, protocol framing characters may be added to the string before transmission."""
    raise HTTPException(status_code=400, detail="Method not implemented.")


@router.put(
    "/api/v1/{device_type}/{device_number}/commandstring",
    responses={
        200: {"model": StringResponse, "description": "Transaction complete or exception"},
        400: {"model": str, "description": "Method or parameter value error, check error message"},
        500: {"model": str, "description": "Server internal error, check error message"},
    },
    tags=["ASCOM Methods Common To All Devices"],
    summary="Transmits an arbitrary string to the device and returns a string value from the device.",
    response_model_by_alias=True,
)
async def device_type_device_number_commandstring_put(
    device_type: str = Path(..., description="One of the recognised ASCOM device types e.g. telescope (must be lower case)", regex=r"^[a-z]*$"),
    device_number: int = Path(..., description="Zero based device number as set on the server (0 to 4294967295)", ge=0, le=4294967295),
    device_type_device_number_commandblind_put_request: DeviceTypeDeviceNumberCommandblindPutRequest = Body(None, description=""),
) -> StringResponse:
    """Transmits an arbitrary string to the device and waits for a string response. Optionally, protocol framing characters may be added to the string before transmission."""
    raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/{device_type}/{device_number}/connected",
    responses={
        200: {"model": BoolResponse, "description": "Driver response"},
        400: {"model": str, "description": "Method or parameter value error, check error message"},
        500: {"model": str, "description": "Server internal error, check error message"},
    },
    tags=["ASCOM Methods Common To All Devices"],
    summary="Retrieves the connected state of the device",
    response_model_by_alias=True,
)
async def device_type_device_number_connected_get(
    device_type: str = Path(..., description="One of the recognised ASCOM device types e.g. telescope (must be lower case)", regex=r"^[a-z]*$"),
    device_number: int = Path(..., description="Zero based device number as set on the server (0 to 4294967295)", ge=0, le=4294967295),
    ClientID: int = Query(1, description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.", ge=0, le=4294967295),
    ClientTransactionID: int = Query(1234, description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.", ge=0, le=4294967295),
) -> BoolResponse:
    """Retrieves the connected state of the device"""
    raise HTTPException(status_code=400, detail="Method not implemented.")


@router.put(
    "/api/v1/{device_type}/{device_number}/connected",
    responses={
        200: {"model": MethodResponse, "description": "Transaction complete or exception"},
        400: {"model": str, "description": "Method or parameter value error, check error message"},
        500: {"model": str, "description": "Server internal error, check error message"},
    },
    tags=["ASCOM Methods Common To All Devices"],
    summary="Sets the connected state of the device",
    response_model_by_alias=True,
)
async def device_type_device_number_connected_put(
    device_type: str = Path(..., description="One of the recognised ASCOM device types e.g. telescope (must be lower case)", regex=r"^[a-z]*$"),
    device_number: int = Path(..., description="Zero based device number as set on the server (0 to 4294967295)", ge=0, le=4294967295),
    Connected: bool = Form(None, description="Set True to connect to the device hardware, set False to disconnect from the device hardware"),
    ClientID: int = Form(None, description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client."),
    ClientTransactionID: int = Form(None, description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs."),
) -> MethodResponse:
    """Sets the connected state of the device"""
    print(f"Connected: {device_type}/{device_number} {Connected}")

    return MethodResponse(ClientTransactionID=ClientTransactionID, ServerTransactionID=ClientTransactionID, ErrorNumber=0, ErrorMessage="")


@router.get(
    "/api/v1/{device_type}/{device_number}/description",
    responses={
        200: {"model": StringResponse, "description": "Driver response"},
        400: {"model": str, "description": "Method or parameter value error, check error message"},
        500: {"model": str, "description": "Server internal error, check error message"},
    },
    tags=["ASCOM Methods Common To All Devices"],
    summary="Device description",
    response_model_by_alias=True,
)
async def device_type_device_number_description_get(
    device_type: str = Path(..., description="One of the recognised ASCOM device types e.g. telescope (must be lower case)", regex=r"^[a-z]*$"),
    device_number: int = Path(..., description="Zero based device number as set on the server (0 to 4294967295)", ge=0, le=4294967295),
    ClientID: int = Query(1, description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.", ge=0, le=4294967295),
    ClientTransactionID: int = Query(1234, description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.", ge=0, le=4294967295),
) -> StringResponse:
    """The description of the device"""
    if device_type == "camera":
        return StringResponse(Value="camera", ClientTransactionID=ClientTransactionID, ServerTransactionID=ClientTransactionID, ErrorNumber=0, ErrorMessage="")



@router.get(
    "/api/v1/{device_type}/{device_number}/driverinfo",
    responses={
        200: {"model": StringResponse, "description": "Driver response"},
        400: {"model": str, "description": "Method or parameter value error, check error message"},
        500: {"model": str, "description": "Server internal error, check error message"},
    },
    tags=["ASCOM Methods Common To All Devices"],
    summary="Device driver description",
    response_model_by_alias=True,
)
async def device_type_device_number_driverinfo_get(
    device_type: str = Path(..., description="One of the recognised ASCOM device types e.g. telescope (must be lower case)", regex=r"^[a-z]*$"),
    device_number: int = Path(..., description="Zero based device number as set on the server (0 to 4294967295)", ge=0, le=4294967295),
    ClientID: int = Query(1, description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.", ge=0, le=4294967295),
    ClientTransactionID: int = Query(1234, description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.", ge=0, le=4294967295),
) -> StringResponse:
    """The description of the driver"""
    raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/{device_type}/{device_number}/driverversion",
    responses={
        200: {"model": StringResponse, "description": "Driver response"},
        400: {"model": str, "description": "Method or parameter value error, check error message"},
        500: {"model": str, "description": "Server internal error, check error message"},
    },
    tags=["ASCOM Methods Common To All Devices"],
    summary="Driver Version",
    response_model_by_alias=True,
)
async def device_type_device_number_driverversion_get(
    device_type: str = Path(..., description="One of the recognised ASCOM device types e.g. telescope (must be lower case)", regex=r"^[a-z]*$"),
    device_number: int = Path(..., description="Zero based device number as set on the server (0 to 4294967295)", ge=0, le=4294967295),
    ClientID: int = Query(1, description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.", ge=0, le=4294967295),
    ClientTransactionID: int = Query(1234, description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.", ge=0, le=4294967295),
) -> StringResponse:
    """A string containing only the major and minor version of the driver."""
    if device_type == "camera":
        return StringResponse(Value="1.0", ClientTransactionID=ClientTransactionID, ServerTransactionID=ClientTransactionID, ErrorNumber=0, ErrorMessage="")
    else:
        raise HTTPException(status_code=400, detail="Method not implemented.")


@router.get(
    "/api/v1/{device_type}/{device_number}/interfaceversion",
    responses={
        200: {"model": IntResponse, "description": "Driver response"},
        400: {"model": str, "description": "Method or parameter value error, check error message"},
        500: {"model": str, "description": "Server internal error, check error message"},
    },
    tags=["ASCOM Methods Common To All Devices"],
    summary="The ASCOM Device interface version number that this device supports.",
    response_model_by_alias=True,
)
async def device_type_device_number_interfaceversion_get(
    device_type: str = Path(..., description="One of the recognised ASCOM device types e.g. telescope (must be lower case)", regex=r"^[a-z]*$"),
    device_number: int = Path(..., description="Zero based device number as set on the server (0 to 4294967295)", ge=0, le=4294967295),
    ClientID: int = Query(1, description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.", ge=0, le=4294967295),
    ClientTransactionID: int = Query(1234, description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.", ge=0, le=4294967295),
) -> IntResponse:
    """This method returns the version of the ASCOM device interface contract to which this device complies. Only one interface version is current at a moment in time and all new devices should be built to the latest interface version. Applications can choose which device interface versions they support and it is in their interest to support  previous versions as well as the current version to ensure thay can use the largest number of devices."""
    print("interface", device_type, device_number)
    if device_type == "camera":
        return IntResponse(Value=3, ClientTransactionID=ClientTransactionID, ServerTransactionID=ClientTransactionID, ErrorNumber=0, ErrorMessage="")


@router.get(
    "/api/v1/{device_type}/{device_number}/name",
    responses={
        200: {"model": StringResponse, "description": "Driver response"},
        400: {"model": str, "description": "Method or parameter value error, check error message"},
        500: {"model": str, "description": "Server internal error, check error message"},
    },
    tags=["ASCOM Methods Common To All Devices"],
    summary="Device name",
    response_model_by_alias=True,
)
async def device_type_device_number_name_get(
    device_type: str = Path(..., description="One of the recognised ASCOM device types e.g. telescope (must be lower case)", regex=r"^[a-z]*$"),
    device_number: int = Path(..., description="Zero based device number as set on the server (0 to 4294967295)", ge=0, le=4294967295),
    ClientID: int = Query(1, description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.", ge=0, le=4294967295),
    ClientTransactionID: int = Query(1234, description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.", ge=0, le=4294967295),
) -> StringResponse:
    """The name of the device"""
    if device_type == "camera":
        return StringResponse(Value="gcamera", ClientTransactionID=ClientTransactionID, ServerTransactionID=ClientTransactionID, ErrorNumber=0, ErrorMessage="")



@router.get(
    "/api/v1/{device_type}/{device_number}/supportedactions",
    responses={
        200: {"model": StringArrayResponse, "description": "Driver response"},
        400: {"model": str, "description": "Method or parameter value error, check error message"},
        500: {"model": str, "description": "Server internal error, check error message"},
    },
    tags=["ASCOM Methods Common To All Devices"],
    summary="Returns the list of action names supported by this driver.",
    response_model_by_alias=True,
)
async def device_type_device_number_supportedactions_get(
    device_type: str = Path(..., description="One of the recognised ASCOM device types e.g. telescope (must be lower case)", regex=r"^[a-z]*$"),
    device_number: int = Path(..., description="Zero based device number as set on the server (0 to 4294967295)", ge=0, le=4294967295),
    ClientID: int = Query(1, description="Client's unique ID. (0 to 4294967295). The client should choose a value at start-up, e.g. a random value between 0 and 65535, and send this value on every transaction to help associate entries in device logs with this particular client.", ge=0, le=4294967295),
    ClientTransactionID: int = Query(1234, description="Client's transaction ID. (0 to 4294967295). The client should start this count at 1 and increment by one on each successive transaction. This will aid associating entries in device logs with corresponding entries in client side logs.", ge=0, le=4294967295),
) -> StringArrayResponse:
    """Returns the list of action names supported by this driver."""
    raise HTTPException(status_code=400, detail="Method not implemented.")
