
# Gaia Camera Simulator

![Example image](example.jpg)

This repositry creates a basic [ASCOM](https://www.ascom-standards.org/) Alpaca Camera driver that simulates the night sky using the [Gaia DR2 catalog](https://www.cosmos.esa.int/web/gaia/dr2) and its photometry.  It requires an Alpaca based telescope driver to be running and connected to the simulator.

When an image is requested, the simulator will:
1. request the telescope's current (ra, dec) coordinates.
2. Use [cabaret](https://github.com/ppp-one/cabaret) to
   1. call Gaia to get the stars in the field of view.
   2. calculate the flux of each star in the field of view.
   3. add noise to the image based on the photon flux from stars, sky background, and camera noise.
3. return the image to the client.

The simulator also has the option of adding scattered sunlight to the image for simulating sky flats. Similarly, one can add (ra, dec) coordinate errors to simulate pointing errors and poor tracking.

## Installation
```shell
git clone https://github.com/ppp-one/gaia-camera-simulator.git
cd gaia-camera-simulator
pip install .
```

## Usage

First start the alpaca telescope simulator -- see [simulators](https://github.com/ASCOMInitiative/ASCOM.Alpaca.Simulators) or [ASCOM Remote](https://github.com/ASCOMInitiative/ASCOMRemote).

Next, adapt the `global_config` such that it corresponds to the IP addresses, device numbers and ports you have configured in the simulator.
```python
from gaia_camera_simulator import ConfigManager

config_manager = ConfigManager()
print(config_manager.global_config)

# Apply changes, e.g. config_manager.global_config["CAMERA_PORT"] = "localhost",
# save using config_manager.save_global_config().
# Alternatively, edit the yaml file directly, located at config_manager.GLOBAL_CONFIG_PATH
```

Now, you can start the camera simulator in a subprocess, e.g. using
```python
from gaia_camera_simulator import GaiaCameraSimulator

gcs = GaiaCameraSimulator()

gcs.start()  # start the simulation subprocess
print(f"Simulation is running {gcs.is_running()}")

telescope, focuser, camera = gcs.initialize_alpaca_devices_from_config(connect_devices=True)
print(f"RA: {telescope.RightAscension}, DEC: {telescope.Declination}")


# Get a simulated image
camera.StartExposure(Duration=1.0, Light=True)
while not camera.ImageReady:
    print("Waiting for image to be ready...")
    time.sleep(0.5)

img = camera.ImageArray  # get image

gcs.stop()  # stop the simulation subprocess

# or, as a context manager
with GaiaCameraSimulator(observatory_config_name="andor_visible") as simulator:
    print(f"Simulation is running {simulator.is_running()}")
    # perform some tests ...
```

Alternatively, you can leave the simulator running for an extended time period as follows:
```shell
python src/gaia_camera_simulator/run.py --host 0.0.0.0 --port 8080 --log-level warning --config_name andor_visible
```

## Configuration

The `ConfigManager` class allows you to create and manage custom observatory configurations for the Gaia camera simulator provided by cabaret.

To add a custom configuration, simply create a `.yaml` file and place it in the directory specified by `config_manager.global_config_dir`. If you wish to change the default directory for storing these configuration files, you can easily do so by following the example below:
```python
from gaia_camera_simulator import ConfigManager

config_manager = ConfigManager()
config_manager.observatory_config_dir = '/path/to/custom/observatory_config_dir/'
config_manager.save_global_config()
print(config_manager)
```

You can store a variation of one of the example configurations as follows:
```python
from gaia_camera_simulator import ConfigManager

config_manager = ConfigManager("default")
config_manager.observatory_config["TRACKING_ERROR"] = True
config_manager.save_observatory_config("default_with_tracking_error")
config_manager.load_observatory_config(config_name="default_with_tracking_error")
print(config_manager.observatory_config["TRACKING_ERROR"])
```
