# nuclides_data.py

nuclides = {
    "130Ba": {
        "name": "Barium 130",
        "half_life": 1.2e21,
        "ideal_daughter": None
    },
    "209Bi": {
        "name": "Bismuth 209",
        "half_life": 2.01e19,
        "ideal_daughter": None
    },
    "113Cd": {
        "name": "Cadmium 113",
        "half_life": 7.70e15,
        "ideal_daughter": None
    },
    "116Cd": {
        "name": "Cadmium 116",
        "half_life": 3.1e19,
        "ideal_daughter": None
    },
    "48Ca": {
        "name": "Calcium 48",
        "half_life": 2.30e19,
        "ideal_daughter": None
    },
    "151Eu": {
        "name": "Europium 151",
        "half_life": 5.00e18,
        "ideal_daughter": None
    },
    "76Ge": {
        "name": "Germanium 76",
        "half_life": 1.8e21,
        "ideal_daughter": None
    },
    "174Hf": {
        "name": "Hafnium 174",
        "half_life": 2.00e15,
        "ideal_daughter": None
    },
    "115In": {
        "name": "Indium 115",
        "half_life": 4.40e14,
        "ideal_daughter": None
    },
    "78K": {
        "name": "Krypton 78",
        "half_life": 9.2e21,
        "ideal_daughter": None
    },
    "100Mo": {
        "name": "Molybdenum 100",
        "half_life": 7.80e18,
        "ideal_daughter": None
    },
    "114Nd": {
        "name": "Neodymium 144",
        "half_life": 2.29e15,
        "ideal_daughter": None
    },
    "150Nd": {
        "name": "Neodymium 150",
        "half_life": 7.90e18,
        "ideal_daughter": None
    },
    "186Os": {
        "name": "Osmium 186",
        "half_life": 2.00e15,
        "ideal_daughter": None
    },
    "148Sm": {
        "name": "Samarium 148",
        "half_life": 7.00e15,
        "ideal_daughter": None
    },
    "82Se": {
        "name": "Selenium 82",
        "half_life": 1.1e20,
        "ideal_daughter": None
    },
    "128Te": {
        "name": "Tellurium 128",
        "half_life": 2.2e24,
        "ideal_daughter": None
    },
    "130Te": {
        "name": "Tellurium 130",
        "half_life": 8.80e18,
        "ideal_daughter": None
    },
    "180W": {
        "name": "Tungsten 180",
        "half_life": 1.80e18,
        "ideal_daughter": None
    },
    "50V": {
        "name": "Vanadium 50",
        "half_life": 1.40e18,
        "ideal_daughter": None
    },
    "124Xe": {
        "name": "Xenon 124",
        "half_life": 1.8e22,
        "ideal_daughter": None
    },
    "136Xe": {
        "name": "Xenon 136",
        "half_life": 2.38e21,
        "ideal_daughter": None
    },
    "96Zr": {
        "name": "Zirconium 96",
        "half_life": 2.00e19,
        "ideal_daughter": None
    },
    "238U": {
        "name": "Uranium 238",
        "half_life": 4.468e9,  # years
        "ideal_daughter": "234Th",
        "daughter_half_life": 6.6e-2,  # 0.066 years ~ 24 days
        "suggested_N0": 1_000_000,
        "suggested_steps": 5000,
        "suggested_time_multiplier": 5
    },
    "235U": {
        "name": "Uranium 235",
        "half_life": 7.04e8,  # years
        "ideal_daughter": "231Pa",
        "daughter_half_life": 0.0013,  # years (~ 11 hours)
        "suggested_N0": 1_000_000,
        "suggested_steps": 5000,
        "suggested_time_multiplier": 5
    },
    "232Th": {
        "name": "Thorium 232",
        "half_life": 1.405e10,  # years
        "ideal_daughter": "228Ra",
        "daughter_half_life": 5.75,   # years
        "suggested_N0": 1_000_000,
        "suggested_steps": 5000,
        "suggested_time_multiplier": 5
    }
}

