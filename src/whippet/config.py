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

    molecules: List[str] = Field(None, description="The list of molecules")

    num_particles: int = Field(0, description="The total number of particles")

    geometry: List[str] = Field(None, description="The geometry to use (plane, pillar)")

    num_images: List[int] = Field(None, description="The number of images to simulate")

    pixel_size: float = Field(1, description="The pixel size to use (A)")

    final_binning: int = Field(
        8, description="The binning for the output reconstruction"
    )


def default() -> Config:
    """
    Return:
        obj: the default configuration

    """
    return Config()


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
