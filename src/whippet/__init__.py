import logging
import os
import parakeet.config
import parakeet.command_line.analyse
import parakeet.command_line.sample
import parakeet.command_line.simulate
from whippet.config import Config


__all__ = ["simulate"]


logger = logging.getLogger(__name__)


def config_initialise(
    config, molecules, num_particles, geometry, num_images, pixel_size
):
    shape = {"plane": "cuboid", "pillar": "cylinder"}[geometry]

    num_molecules = num_particles

    num_images = num_images
    dose = 140 / num_images

    if shape == "cuboid":
        start_angle = -60
        step_angle = 120 / num_images
    else:
        start_angle = -90
        step_angle = 180 / num_images

    num_particles_per_molecule = num_particles // len(molecules)
    pdb = [
        {"id": pdb_id, "instances": num_particles_per_molecule} for pdb_id in molecules
    ]

    # Write the config file
    config_obj = parakeet.config.Config(
        **{
            "microscope": {
                "beam": {
                    "electrons_per_angstrom": dose,
                    "energy": 300,
                    "illumination_semiangle": 0.02,
                    "acceleration_voltage_spread": 8.0e-07,
                    "energy_spread": 2.66e-06,
                },
                "detector": {
                    "nx": 4000,
                    "ny": 4000,
                    "pixel_size": pixel_size,
                    "origin": [3000, 3000],
                },
                "lens": {
                    "c_10": -30000,
                    "c_30": 2.7,
                    "c_c": 2.7,
                    "current_spread": 3.3e-07,
                },
            },
            "sample": {
                "box": [12000, 12000, 12000],
                "centre": [6000, 6000, 6000],
                "shape": {
                    "cuboid": {
                        "length_x": 6000,
                        "length_y": 6000,
                        "length_z": 1500,
                    },
                    "cylinder": {
                        "length": 8000,
                        "radius": 1500,
                    },
                    "type": shape,
                },
                "molecules": {"pdb": pdb},
            },
            "scan": {
                "axis": [0, 1, 0],
                "exposure_time": 1,
                "mode": "tilt_series",
                "num_images": num_images,
                "start_angle": start_angle,
                "start_pos": 0,
                "step_angle": step_angle,
                "step_pos": "auto",
            },
            "simulation": {
                "division_thickness": 10000,
                "ice": True,
                "inelastic_model": None,
                "mp_loss_width": None,
                "margin": 100,
                "padding": 100,
                "radiation_damage_model": False,
                "sensitivity_coefficient": 0.014,
                "slice_thickness": 5.0,
            },
        }
    )

    parakeet.config.save(config_obj, config)


def simulate_and_reconstruct_single(
    molecules, num_particles, geometry, num_images, pixel_size, final_binning
):

    # The data directory
    data_directory = os.path.join(geometry, "%d" % num_images)
    if not os.path.exists(data_directory):
        os.makedirs(data_directory)

    # Set filenames
    config = os.path.join(data_directory, "config.yaml")
    sample = os.path.join(data_directory, "sample.h5")
    exit_wave = os.path.join(data_directory, "exit_wave.mrc")
    optics = os.path.join(data_directory, "optics.mrc")
    image = os.path.join(data_directory, "image.mrc")
    config_rebinned = os.path.join(data_directory, "config.yaml")
    rec = os.path.join(data_directory, "rec.mrc")

    # Setup the config file
    if not os.path.exists(config):
        config_initialise(
            config, molecules, num_particles, geometry, num_images, pixel_size
        )

    # Generate the sample
    if not os.path.exists(sample):
        parakeet.command_line.sample.new(["-c", config, "-s", sample])
        parakeet.command_line.sample.add_molecules(["-c", config, "-s", sample])

    # Simulate the exit wave
    if not os.path.exists(exit_wave):
        parakeet.command_line.simulate.exit_wave(
            ["-c", config, "-s", sample, "-e", exit_wave]
        )

    # Simulate optics
    if not os.path.exists(optics):
        parakeet.command_line.simulate.optics(
            ["-c", config, "-e", exit_wave, "-o", optics]
        )

    # Simulate image
    if not os.path.exists(image):
        parakeet.command_line.simulate.image(["-c", config, "-o", optics, "-i", image])

    # Setup the config file
    if not os.path.exists(config):
        config_initialise(
            config,
            molecules,
            num_particles,
            geometry,
            num_images,
            pixel_size * final_binning,
        )

    # Do the reconstruction
    if not os.path.exists(rec):
        parakeet.command_line.analyse.reconstruct(
            ["-c", config, "-i", image, "-r", rec]
        )


def simulate_and_reconstruct(config: Config):
    """
    Simulate and reconstruct

    """
    # Loop through all the parameters and do the simulation
    for geometry in config.geometry:
        for num_images in config.num_images:
            simulate_and_reconstruct_single(
                config.molecules,
                config.num_particles,
                geometry,
                num_images,
                config.pixel_size,
                config.final_binning,
            )
