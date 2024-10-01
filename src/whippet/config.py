#
# This code is distributed under the GPLv3 license, a copy of
# which is included in the root directory of this package.
#
import logging
import yaml

from pydantic import BaseModel as PydanticBaseModel
from pydantic import Field
from typing import List
from typing import Union

# Get the logger
logger = logging.getLogger(__name__)


class BaseModel(PydanticBaseModel):
    """
    Create a custom base to define desired behaviour

    """

    class Config:
        # Ensure that enums use string values
        use_enum_values = True

        # Don't allow extra fields
        extra = "forbid"


class Config(BaseModel):
    """
    The whippet configuration parameters

    """

    pixel_size: float = Field(1.9, description="The pixel size to use (A)")

    final_binning: int = Field(
        8, description="The binning for the output reconstruction"
    )

    pdb: List[str] = Field(None, description="The PDB filenames")

    tile_size: float = Field(None, description="The tiling size (A)")

    tile_angle: float = Field(0.0, description="The tiling angle (degrees)")

    all_orientations: bool = Field(
        False, description="Simulate lamellae with all orientations"
    )

    thickness: float = Field(1000, description="The thickness of the lamella (A)")

    start_angle: float = Field(51, description="The starting tilt angle (deg)")

    step_angle: float = Field(3, description="The tilt angle step (deg)")

    num_images: int = Field(1, description="The number of images")


def default() -> Config:
    """
    Return:
        obj: the default configuration

    """
    return Config()


def example() -> Config:
    """
    Return:
        The configuration object

    """
    return Config(
        **{
            "pdb": ["my_file.pdb"],
            "tile_size": 200,
            "tile_angle": 5,
            "all_orientations": True,
            "thickness": 1000,
            "start_angle": -51,
            "step_angle": 3,
            "num_images": 35,
            "pixel_size": 1.9,
            "final_binning": 8,
        }
    )


def save(config: Config, filename: str = "config.yaml", **kwargs):
    """
    Save the configuration file

    Args:
        config (str): The configuration object
        filename (str): The configuration filename

    """

    # Get the dictionary
    d = config.dict(**kwargs)

    # Write the output file
    with open(filename, "w") as outfile:
        yaml.safe_dump(d, outfile)


def load(config: Union[str, dict] = None) -> Config:
    """
    Load the configuration from the various inputs

    Args:
        config (str): The config filename or config dictionary

    Returns:
        dict: The configuration dictionary

    """

    # If the yaml configuration is set then merge the configuration
    if config:
        if isinstance(config, str):
            with open(config) as infile:
                config_file = yaml.safe_load(infile)
        else:
            config_file = config
    else:
        config_file = {}

    # Get the configuration
    return Config(**config_file)
