import logging
import mrcfile
import numpy as np
import os
import parakeet.config
import parakeet.command_line.analyse
import parakeet.command_line.sample
import parakeet.command_line.simulate
from whippet.config import Config
from math import sqrt, ceil


__all__ = ["simulate"]


logger = logging.getLogger(__name__)


def config_initialise(
    config, molecules, num_particles, geometry, num_images, pixel_size, final_binning=1
):
    shape = {"plane": "cuboid", "pillar": "cylinder"}[geometry]

    num_molecules = num_particles

    num_images = num_images
    dose = 140 / num_images

    if shape == "cuboid":
        start_angle = -60
        step_angle = 120 / num_images
        thickness = 1500
        d2 = 4000 * pixel_size / 2.0
        l2 = thickness / 2.0
        o2 = 12000 - 2000 * pixel_size
        margin = [
            max(0, d2 - d2 / sqrt(2)) + o2 + 100,
            o2 + 100,
            max(0, l2 - d2 / sqrt(2)) + 100,
        ]
    else:
        start_angle = -90
        step_angle = 180 / num_images
        thickness = 1500 * 2
        d2 = 4000 * pixel_size / 2.0
        l2 = thickness / 2.0
        o2 = 12000 - 2000 * pixel_size
        margin = [max(0, l2 - d2), o2 + 100, max(0, l2 - d2)]

    image_size = 4000 // final_binning
    pixel_size = pixel_size * final_binning

    origin = [12000 - 2000 * pixel_size, 12000 - 2000 * pixel_size]

    num_particles_per_molecule = int(ceil(num_particles / len(molecules)))
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
                    "nx": image_size,
                    "ny": image_size,
                    "pixel_size": pixel_size,
                    "origin": origin,
                },
                "lens": {
                    "c_10": -30000,
                    "c_30": 2.7,
                    "c_c": 2.7,
                    "current_spread": 3.3e-07,
                },
            },
            "sample": {
                "box": [24000, 24000, 24000],
                "centre": [12000, 12000, 12000],
                "shape": {
                    "cuboid": {
                        "length_x": 24000,
                        "length_y": 24000,
                        "length_z": 1500,
                    },
                    "cylinder": {
                        "length": 24000,
                        "radius": 1500,
                    },
                    "type": shape,
                    "margin": margin,
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


def write_coordinates(sample_filename, rec_filename, coordinates_filename):

    with open(coordinates_filename, "w") as outfile:

        # Get the size of the volume
        tomo_file = mrcfile.mmap(rec_filename)
        tomogram = tomo_file.data
        voxel_size = np.array(
            (
                tomo_file.voxel_size["x"],
                tomo_file.voxel_size["y"],
                tomo_file.voxel_size["z"],
            )
        )
        assert voxel_size[0] == voxel_size[1]
        assert voxel_size[0] == voxel_size[2]
        size = np.array(tomogram.shape)[[2, 0, 1]] * voxel_size

        # Open sample and get the centre
        sample = parakeet.sample.Sample(sample_filename)
        centre = np.array(sample.centre)

        # Write out the coordinates in (X, Y, Z). The axes of the reconstruction are (Y, Z, X)
        for molecule, (_, positions, orientations) in sample.iter_molecules():
            for position, orientation in zip(positions, orientations):
                position = position - (centre - size / 2.0)
                position[2] = size[2] - position[2]
                outfile.write(
                    "%s, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n"
                    % (
                        molecule,
                        position[0],
                        position[1],
                        position[2],
                        orientation[0],
                        orientation[1],
                        orientation[2],
                    )
                )


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
    config_rebinned = os.path.join(data_directory, "config_rebinned.yaml")
    image_rebinned = os.path.join(data_directory, "image_rebinned.mrc")
    rec = os.path.join(data_directory, "rec.mrc")
    coordinates = os.path.join(data_directory, "coords.csv")

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
    if not os.path.exists(config_rebinned):
        config_initialise(
            config_rebinned,
            molecules,
            num_particles,
            geometry,
            num_images,
            pixel_size,
            final_binning,
        )

    # Rebin the image
    if not os.path.exists(image_rebinned):
        parakeet.command_line.export(
            [image, "-o", image_rebinned, "--rebin=%d" % final_binning]
        )

    # Do the reconstruction
    if not os.path.exists(rec):
        parakeet.command_line.analyse.reconstruct(
            ["-c", config_rebinned, "-i", image_rebinned, "-r", rec]
        )

    # Extract the coordinates
    if not os.path.exists(coordinates):
        write_coordinates(sample, rec, coordinates)


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
