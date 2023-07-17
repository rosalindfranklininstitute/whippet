from whippet.config import Config
from whippet import simulate_and_reconstruct


def test_whippet():

    # Create the configuration
    config = {
        "molecules": ["4v1w"],
        "num_particles": 100,
        "geometry": ["plane"],
        "num_images": [40],
        "pixel_size": 1,
        "final_binning": 8,
    }

    # Simulate and reconstruct
    simulate_and_reconstruct(Config(**config))
