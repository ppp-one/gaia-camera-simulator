import argparse
import logging
from typing import Tuple

import uvicorn

from gaia_camera_simulator.config_manager import ConfigManager

logging.basicConfig(level=logging.INFO)


def parse_arguments() -> Tuple[argparse.Namespace, dict]:
    def pair_up_args(arg_list):
        return {arg_list[i]: arg_list[i + 1] for i in range(0, len(arg_list), 2)}

    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "-c",
        "--config_name",
        type=str,
        help="Name of the configuration file to use. Defaults to 'default'.",
        default="default",
    )
    argparser.add_argument(
        "-cp",
        "--config_path",
        type=str,
        help="Path to the configuration file to use. Defaults to None.",
        default="",
    )

    args, unknown_args_list = argparser.parse_known_args()
    return args, pair_up_args(unknown_args_list)


if __name__ == "__main__":
    args, uvicorn_args = parse_arguments()
    config_name = str(args.config_name)
    config_file_path = str(args.config_path)
    if not config_file_path or config_file_path == "None":
        config_file_path = None
    
    ConfigManager(
        observatory_config_name=config_name,
        observatory_config_file_path=config_file_path,
    )

    default_args = {
        "host": "0.0.0.0",
        "port": 8080,
        "log_level": "warning",
    }
    uvicorn_args = {
        key: value
        for key, value in uvicorn_args.items()
        if key in uvicorn.run.__code__.co_varnames
    }

    uvicorn_args = default_args | uvicorn_args

    uvicorn.run("main:app", **uvicorn_args)
