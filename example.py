import time

import matplotlib.pyplot as plt
import numpy as np
from alpaca.camera import ImageArrayElementTypes
from astropy.visualization import ZScaleInterval

from gaia_camera_simulator import GaiaCameraSimulator

# Use the ConfigManager to set save and load custom configurations
gcs = GaiaCameraSimulator(
    observatory_config_name="default", host="localhost", port=8080
)
config = gcs.config

# Start simulation subprocess
gcs.start()

print("Running example with config:")
print(config)

telescope, focuser, camera = gcs.initialize_alpaca_devices_from_config(connect_devices=True)

# Move focuser to maximal focus position 
focuser.Move(config.get("FOCUSER_POSITION", 10_000))
while focuser.IsMoving:
    print("Waiting for focuser to stop moving...")
    time.sleep(0.5)

# move to target (M2)
# turn tracking on
telescope.Tracking = True
ra = 323.36152
dec = -0.82325

telescope.SlewToCoordinatesAsync(RightAscension=ra * 24 / 360, Declination=dec)
while telescope.Slewing:
    print("Waiting for slew to complete...")
    time.sleep(0.5)

# start exposure
camera.StartExposure(1.0, True)  # change to False for dark frame

# wait for exposure to complete/image ready
while not camera.ImageReady:
    print("Waiting for image to be ready...")
    time.sleep(0.5)

# get image
img = camera.ImageArray

# get image info such displayed correctly
imginfo = camera.ImageArrayInfo
if imginfo.ImageElementType == ImageArrayElementTypes.Int32:
    if camera.MaxADU <= 65535:
        imgDataType = np.uint16  # Required for BZERO & BSCALE to be written
    else:
        imgDataType = np.int32
elif imginfo.ImageElementType == ImageArrayElementTypes.Double:
    imgDataType = np.float64

# Make a numpy array of he correct shape for astropy.io.fits (if needed)
if imginfo.Rank == 2:
    nda = np.array(img, dtype=imgDataType).transpose()
else:
    nda = np.array(img, dtype=imgDataType).transpose(2, 1, 0)

# display the image
zscale = ZScaleInterval(contrast=0.05, krej=2.5)
vmin, vmax = zscale.get_limits(values=nda)

fig = plt.figure()
plt.imshow(
    nda,
    vmin=vmin,
    vmax=vmax,
    origin="lower",
    cmap=plt.get_cmap("gray"),
    interpolation="none",
)
# plt.axis('off')
# fig.savefig('image.jpg', dpi=300, bbox_inches='tight', pad_inches=0)
plt.show()

# disconnect from telescope
telescope.Connected = False

# disconnect from camera
camera.Connected = False

# Stop simulation subprocess
gcs.stop()