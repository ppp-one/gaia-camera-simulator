import logging
import subprocess
import time
from pathlib import Path

import httpx

from gaia_camera_simulator.config_manager import ConfigManager

SCRIPT_PATH = Path(__file__).parent / "run.py"


class GaiaCameraSimulator:
    """A class to simulate a Gaia camera.


    Example
    -------
    >>> from gaia_camera_simulator import GaiaCameraSimulator
    >>> gcs = GaiaCameraSimulator(observatory_config_name="deault")
    >>> gcs.start()
    >>> gcs.is_running()
    >>> gcs.process.communicate()
    >>> gcs.stop()

    >>> with GaiaCameraSimulator(observatory_config_name="my_config") as simulator:
    >>>     simulator.start()
    >>>     if simulator.is_running():
    >>>     print("Simulation is running.")
    """

    _SCRIPT_PATH = SCRIPT_PATH

    def __init__(
        self,
        observatory_config_name: str = "default",
        observatory_config_file_path: Path | None = None,
        host: str | None= None,
        port: int | None= None,
        log_level: str | None = None,
    ):
        self.observatory_config_name = observatory_config_name
        self.observatory_config_file_path = observatory_config_file_path
        config = self.config

        self.host =  host if host is not None else config.get("CAMERA_HOST", "localhost")
        self.port = port if port is not None else config.get("CAMERA_PORT", "8080")
        self.log_level = log_level if log_level is not None else config.get("LOG_LEVEL", "warning")
        self.process = None

    @property
    def config(self):
        config = ConfigManager(
            observatory_config_name=self.observatory_config_name,
            observatory_config_file_path=self.observatory_config_file_path,
        ).combine_configs()
        return config

    @property
    def stdout(self):
        if self.process is None:
            return None
        return self.process.stdout

    def start(self):
        """Start the simulation process."""
        if self.process is not None:
            logging.info("Simulation is already running.")
            return

        try:
            self.process = subprocess.Popen(
                [
                    "python",
                    str(self._SCRIPT_PATH),
                    "--host",
                    str(self.host),
                    "--port",
                    str(self.port),
                    "--log_level",
                    str(self.log_level),
                    "--config_name",
                    str(self.observatory_config_name),
                    "--config_path",
                    str(self.observatory_config_file_path)
                    if self.observatory_config_file_path
                    else "",
                ],
                stdout=subprocess.PIPE,
            )
            # Check if the server is running
            logging.info("Waiting for the server to start.")
            self.wait_for_server()

            logging.info("Simulation started.")
        except Exception as e:
            logging.error(f"Failed to start simulation: {e}")

    def wait_for_server(self, timeout=30):
        start_time = time.time()
        while True:
            try:
                response = httpx.get(f"http://{self.host}:{self.port}/")
                if response.headers.get("server") == "uvicorn":
                    logging.info("Server started successfully.")
                    return
            except httpx.RequestError as e:
                logging.debug(f"Request error: {e}. Waiting for the server to start.")
            except Exception as e:
                logging.error(
                    f"Unexpected error: {e}. Waiting for the server to start."
                )

            if time.time() - start_time > timeout:
                logging.error("Timeout while waiting for the server to start.")
                self.process.kill()  # Ensure the process is terminated on timeout
                raise RuntimeError("Server did not start in the expected time.")

            time.sleep(0.05)

    def stop(self):
        """Stop the simulation process."""
        if self.process is None:
            logging.warning("No simulation is running.")
            return

        try:
            self.process.terminate()
            self.process.wait()
            logging.info("Simulation stopped.")
        except Exception as e:
            logging.error(f"Failed to stop simulation: {e}")
        finally:
            self.process = None

    def is_running(self) -> bool:
        """Check if the simulation is currently running."""
        return self.process is not None and self.process.poll() is None

    def cleanup(self):
        """Clean up resources."""
        if self.is_running():
            self.stop()
            
    def initialize_alpaca_devices_from_config(self, connect_devices=False):
        """Initialize the Alpaca devices using the configuration."""
        from alpaca.camera import Camera
        from alpaca.focuser import Focuser
        from alpaca.telescope import Telescope
        
        telescope = Telescope(self.config["TELESCOPE_IP"], self.config["TELESCOPE_DEVICE_NUMBER"])
        focuser = Focuser(self.config["FOCUSER_IP"], self.config["FOCUSER_DEVICE_NUMBER"])
        camera = Camera(f"{self.host}:{self.port}", 0)
        
        if connect_devices:
            telescope.Connected = True
            focuser.Connected = True
            camera.Connected = True
        
        return telescope, focuser, camera

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()

    def __del__(self):
        self.cleanup()

    def __repr__(self):
        return (
            f"GaiaCameraSimulator("
            f"observatory_config_name={self.observatory_config_name}, "
            f"observatory_config_file_path={self.observatory_config_file_path}, "
            f"host={self.host}, port={self.port}, log_level={self.log_level},)"
        )