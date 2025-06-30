from pathlib import Path
from shutil import copyfile

import yaml

DEFAULT_OBSERVATORY_CONFIG_DIR = Path(__file__).parent / "observatory_configs"


class ConfigManager:
    """A class to load and save configuration files for the Gaia camera simulator.

    This is a singleton class that loads the global configuration file and the observatory
    configuration file. The global configuration file is located in the same directory as this
    module and is named `global_config.yaml`. The observatory configuration files are located in
    the `observatory_configs` directory in the same directory as this module.

    Parameters
    ----------
    observatory_config_file_path : Path, optional
        Path to the observatory configuration file, by default None
    config_name : str, optional
        Name of the configuration file, by default "default"

    Examples
    --------
    >>> from gaia_camera_simulator import ConfigManager
    >>> config_manager = ConfigManager("default")
    >>> print(config_manager)
    """

    _instance = None
    DEFAULT_OBSERVATORY_CONFIG_DIR = DEFAULT_OBSERVATORY_CONFIG_DIR
    GLOBAL_CONFIG_PATH = Path(__file__).parent / "global_config.yaml"
    observatory_config = {}
    global_config = {}

    def __new__(
        cls,
        observatory_config_name="default",
        observatory_config_file_path=None,
    ):
        if cls._instance is None:
            cls._instance = super(ConfigManager, cls).__new__(cls)
            cls._instance.load_global_config()
            cls._instance.load_observatory_config(
                observatory_config_name,
                observatory_config_file_path,
            )
        return cls._instance

    def combine_configs(self):
        return self.global_config | self.observatory_config

    def __getitem__(self, key):
        if key in self.global_config:
            return self.global_config[key]
        elif key in self.observatory_config:
            return self.observatory_config[key]
        else:
            raise KeyError(f"Key {key} not found in global or observatory config.")

    def __setitem__(self, key, value):
        if key in self.global_config:
            self.global_config[key] = value
        else:
            self.observatory_config[key] = value

    @property
    def observatory_config_dir(self):
        if not self.global_config:
            self.load_global_config()

        if "OBSERVATORY_CONFIG_DIR" not in self.global_config:
            self.observatory_config_dir = str(self.DEFAULT_OBSERVATORY_CONFIG_DIR)

        return Path(self.global_config["OBSERVATORY_CONFIG_DIR"])

    @observatory_config_dir.setter
    def observatory_config_dir(self, value):
        self.global_config["OBSERVATORY_CONFIG_DIR"] = str(value)

    @property
    def observatory_config_name(self):
        return self.observatory_config.get("config_name", "default")

    @property
    def observatory_config_file_path(self):
        return self.observatory_config_dir / f"{self.observatory_config_name}.yaml"

    def load_global_config(self):
        if not self.GLOBAL_CONFIG_PATH.exists():
            self.create_global_config()
            return

        with open(self.GLOBAL_CONFIG_PATH, "r") as stream:
            self.global_config = yaml.safe_load(stream)

    def load_observatory_config(
        self, config_name: str | None, config_file_path: Path | None = None
    ):
        if config_name is None and config_file_path is None:
            config_name = "default"

        if config_file_path is None:
            config_file_path = self.observatory_config_dir / f"{config_name}.yaml"

        if not config_file_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_file_path}.")

        if config_name is None:
            config_name = config_file_path.stem

        with open(config_file_path, "r") as stream:
            observatory_config = yaml.safe_load(stream)
            observatory_config["config_name"] = config_name
            self.observatory_config = observatory_config

    def save_observatory_config(self, config_name, config=None):
        if config is None:
            config = self.observatory_config

        config_file_path = self.observatory_config_dir / f"{config_name}.yaml"

        with open(config_file_path, "w") as stream:
            yaml.dump(config, stream)

    def save_global_config(self):
        with open(self.GLOBAL_CONFIG_PATH, "w") as stream:
            yaml.dump(self.global_config, stream)

        # Copy existing
        if self.observatory_config_dir != self.DEFAULT_OBSERVATORY_CONFIG_DIR:
            if not self.observatory_config_dir.exists():
                self.observatory_config_dir.mkdir(parents=True)

            existing_files = [f.name for f in self.observatory_config_dir.iterdir()]

            for file_name in self.DEFAULT_OBSERVATORY_CONFIG_DIR.iterdir():
                if file_name not in existing_files:
                    copyfile(
                        self.DEFAULT_OBSERVATORY_CONFIG_DIR / file_name,
                        self.observatory_config_dir / file_name,
                    )

    def create_global_config(self):
        self.global_config = {
            "TELESCOPE_IP": "localhost",
            "TELESCOPE_DEVICE_NUMBER": 0,
            "FOCUSER_IP": "localhost",
            "FOCUSER_DEVICE_NUMBER": 0,
            "CAMERA_IP": "localhost",
            "CAMERA_PORT": "8080",
            "LOG_LEVEL": "warning",
            "OBSERVATORY_CONFIG_DIR": str(self.DEFAULT_OBSERVATORY_CONFIG_DIR),
        }
        self.save_global_config()

    def reset_global_config(self):
        self.create_global_config()

    def __repr__(self) -> str:
        return (
            "ConfigManager("
            f"observatory_config={self.observatory_config}, global_config={self.global_config})"
        )

    def __str__(self) -> str:
        return (
            "ConfigManager(\n"
            f"  observatory_config={self.observatory_config},\n"
            f"  global_config={self.global_config}\n"
            ")"
        )
