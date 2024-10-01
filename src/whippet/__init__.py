import gemmi
import logging
import mrcfile
import numpy as np
import os
import parakeet.config
import parakeet.command_line.analyse
import parakeet.command_line.sample
import parakeet.command_line.simulate
from scipy.spatial.transform import Rotation as R
from whippet.config import Config
from math import floor


__all__ = ["simulate"]


logger = logging.getLogger(__name__)


def config_initialise(
    config,
    filename,
    orientation,
    tile_size,
    tile_angle,
    thickness,
    start_angle,
    step_angle,
    num_images,
    pixel_size,
    final_binning=1,
):
    def AxB(A, B):
        return R.from_matrix(np.matmul(A.as_matrix(), B.as_matrix())).as_rotvec()

    # Tilt the tiling pattern to avoid aligning exactly along pixels
    rotation = R.from_rotvec((0, 0, np.radians(tile_angle)))

    # Compute the orientation to put each face up
    orientation = {
        # no rotation
        0: np.array((0, 0, 0)),
        #   rotation of 90 degrees around the X-axis
        1: np.array((1, 0, 0)) * np.pi / 2,
        # rotation of -90 degrees around the X-axis
        2: np.array((1, 0, 0)) * -np.pi / 2,
        # rotation of 90 degrees around the Y-axis
        3: np.array((0, 1, 0)) * np.pi / 2,
        # rotation of -90 degrees around the Y-axis
        4: np.array((0, 1, 0)) * -np.pi / 2,
        # rotation of 180 degrees around the X-axis
        5: np.array((1, 0, 0)) * 2 * np.pi / 2,
    }[orientation]

    # Apply the rotation to the orientation
    orientation = AxB(rotation, R.from_rotvec(orientation))

    # Split the dose over the number of images
    dose = 175 / num_images

    # Compute the image size and pixel size
    image_size = 4096 // final_binning
    pixel_size = pixel_size * final_binning

    # Compute the lenght and centre of the sample space
    length_x = 24000
    length_y = 24000
    centre_x = length_x / 2
    centre_y = length_y / 2
    centre_z = length_x / 2
    origin = [
        centre_x - 0.5 * image_size * pixel_size,
        centre_y - 0.5 * image_size * pixel_size,
    ]

    # Do the tiling of the input PDB file
    num_tiles_x = int(floor(length_x / tile_size))
    num_tiles_y = int(floor(length_y / tile_size))
    tiles = []
    for j in range(num_tiles_y):
        for i in range(num_tiles_x):
            x = i * tile_size
            y = j * tile_size
            z = centre_z
            p = np.matmul(rotation.as_matrix(), (x, y, z))

            # Only include particles within field of view
            if (p[1] > centre_y - 0.5 * image_size * pixel_size) and (
                p[1] < centre_y + 0.5 * image_size * pixel_size
            ):
                tiles.append(
                    {
                        "position": tuple(map(float, p)),
                        "orientation": tuple(map(float, orientation)),
                    }
                )

    # Write the config file
    config_obj = parakeet.config.Config(
        **{
            "microscope": {
                "beam": {
                    "electrons_per_angstrom": dose,
                    "energy": 300,                   #
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
                "box": [length_x, length_y, length_x],
                "centre": [centre_x, centre_y, centre_z],
                "shape": {
                    "cuboid": {
                        "length_x": length_x,
                        "length_y": length_y,
                        "length_z": thickness,
                    },
                    "type": "cuboid",
                },
                "molecules": {
                    "local": [
                        {
                            "filename": filename,
                            "instances": tiles,
                        }
                    ],
                },
            },
            "scan": {
                "axis": [0, 1, 0],
                "exposure_time": 1.67,
                "mode": "dose_symmetric",
                "num_images": num_images,
                "start_angle": start_angle,
                "start_pos": 0,
                "step_angle": step_angle,
                "step_pos": "auto",
            },
            "simulation": {
                "division_thickness": 10000,
                "ice": True,                       #
                "inelastic_model": None,
                "mp_loss_width": None,
                "margin": 100,
                "padding": 100,
                "radiation_damage_model": True,    #
                "sensitivity_coefficient": 0.022,   #
                "slice_thickness": 5.0,
            },
        }
    )

    parakeet.config.save(config_obj, config)


def compute_centre_of_mass(pdb_filename):

    structure = gemmi.read_structure(pdb_filename)
    centre_of_mass = {}
    for model in structure:
        global_centre_of_mass = np.array(model.calculate_center_of_mass().tolist())
        for chain in model:
            centre_of_mass[chain.name] = list(
                np.array(chain.calculate_center_of_mass().tolist())
                - global_centre_of_mass
            )

    return centre_of_mass


def write_coordinates(
    sample_filename, rec_filename, coordinates_filename, pdb_filename
):

    # Compute the centre of mass of the particles within the tiles
    centre_of_mass = compute_centre_of_mass(pdb_filename)
    for chain, com in centre_of_mass.items():
        print(
            "Centre of mass of %s = (%.2f, %.2f, %.2f)"
            % (chain, com[0], com[1], com[2])
        )

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
                print(position, centre, size / 2.0)
                position = position - (centre - size / 2.0)
                for chain, com in centre_of_mass.items():

                    # Rotate the relative position to find the position of the
                    # particle within the volume
                    p = np.matmul(R.from_rotvec(orientation).as_matrix(), com)
                    p = position + p
                    p[2] = size[2] - p[2]
                    p = p / voxel_size

                    # Only print out particles in the volume
                    if np.all((p > 0) & (p < size)):
                        outfile.write(
                            "%s, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n"
                            % (
                                chain,
                                p[0],
                                p[1],
                                p[2],
                                orientation[0],
                                orientation[1],
                                orientation[2],
                            )
                        )


def compute_particle_size(pdb_filename):

    structure = gemmi.read_structure(pdb_filename)
    particle_size = {}
    for model in structure:
        for chain in model:
            centre = np.array(chain.calculate_center_of_mass().tolist())
            atoms = []
            for residue in chain:
                for atom in residue:
                    atoms.append(atom.pos.tolist())
            distance = np.sqrt(np.sum((centre - np.array(atoms)) ** 2, axis=1))
            particle_size[chain.name] = np.max(distance)
    return particle_size


def average_particles(rec, coordinates, pdb_filename, average_prefix):

    # Get the reconstruction
    tomo_file = mrcfile.open(rec)
    voxel_size = np.array(
        (
            tomo_file.voxel_size["x"],
            tomo_file.voxel_size["y"],
            tomo_file.voxel_size["z"],
        )
    )
    assert voxel_size[0] == voxel_size[1]
    assert voxel_size[0] == voxel_size[2]

    particle_size = compute_particle_size(pdb_filename)

    average = {}
    count = {}
    for name, size in particle_size.items():
        shape = np.ceil(
            np.array([1 * 1 * 2 * size * np.sqrt(2)] * 3) / voxel_size
        ).astype(int)
        average[name] = np.zeros(shape)
        count[name] = 0

    with open(coordinates) as infile:
        for line in infile.readlines():
            tokens = line.split(",")
            name = tokens[0]
            p = np.array(list(map(float, tokens[1:4])))
            r = np.array(list(map(float, tokens[4:])))

            shape = average[name].shape
            x0 = int(np.floor(p[0] - shape[2] // 2))
            x1 = int(np.floor(x0 + shape[2]))
            y0 = int(np.floor(p[1] - shape[1] // 2))
            y1 = int(np.floor(y0 + shape[1]))
            z0 = int(np.floor(p[2] - shape[0] // 2))
            z1 = int(np.floor(z0 + shape[0]))

            if (
                (x0 >= 0)
                and (y0 >= 0)
                and (z0 >= 0)
                and (x1 <= tomo_file.data.shape[2])
                and (z1 <= tomo_file.data.shape[1])
                and (y1 <= tomo_file.data.shape[0])
            ):
                particle = tomo_file.data[y0:y1, z0:z1, x0:x1]

                average[name] += particle
                count[name] += 1

    for name in average:
        filename = "%s_%s.mrc" % (average_prefix, name)
        if not os.path.exists(filename):
            print("%s: writing average of %d particles" % (name, count[name]))
            data = average[name] / count[name]
            handle = mrcfile.new(filename, overwrite=True)
            handle.set_data(data.astype("float32"))
            del handle


def simulate_and_reconstruct_single(
    filename,
    orientation,
    tile_size,
    tile_angle,
    thickness,
    start_angle,
    step_angle,
    num_images,
    pixel_size,
    final_binning,
):

    # The data directory
    identifier = os.path.basename(filename)
    data_directory = os.path.abspath(
        os.path.join("sim", "%s_face_%d" % (identifier, orientation))
    )
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
    average_prefix = os.path.join(data_directory, "average")

    # Setup the config file
    if not os.path.exists(config):
        config_initialise(
            config,
            filename,
            orientation,
            tile_size,
            tile_angle,
            thickness,
            start_angle,
            step_angle,
            num_images,
            pixel_size,
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
            filename,
            orientation,
            tile_size,
            tile_angle,
            thickness,
            start_angle,
            step_angle,
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
        write_coordinates(sample, rec, coordinates, filename)

    # Average the particles
    average_particles(rec, coordinates, filename, average_prefix)


def simulate_and_reconstruct(config: Config):
    """
    Simulate and reconstruct

    """

    def get_orientations(all_orientations):
        return list(range(6)) if all_orientations else [0]

    # Loop through all the parameters and do the simulation
    for filename in config.pdb:
        for orientation in get_orientations(config.all_orientations):
            simulate_and_reconstruct_single(
                filename,
                orientation,
                config.tile_size,
                config.tile_angle,
                config.thickness,
                config.start_angle,
                config.step_angle,
                config.num_images,
                config.pixel_size,
                config.final_binning,
            )
