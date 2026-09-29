# Half-lives in years, from NUBASE2020:
#   F.G. Kondev et al., Chinese Physics C 45, 030001 (2021).
# 130Ba, 78Kr and 124Xe are listed as stable in NUBASE2020; the values below
# are measured double electron capture half-lives.

nuclides = {
    "130Ba": {
        "name": "Barium 130",
        "half_life": 1.2e+21,  # double electron capture; NUBASE2020: ~1e21 y
        "daughter": None
    },
    "209Bi": {
        "name": "Bismuth 209",
        "half_life": 2.01e+19,  # alpha
        "daughter": None
    },
    "113Cd": {
        "name": "Cadmium 113",
        "half_life": 8.04e+15,  # beta-
        "daughter": None
    },
    "116Cd": {
        "name": "Cadmium 116",
        "half_life": 2.69e+19,  # double beta-
        "daughter": None
    },
    "48Ca": {
        "name": "Calcium 48",
        "half_life": 5.6e+19,  # double beta-
        "daughter": None
    },
    "151Eu": {
        "name": "Europium 151",
        "half_life": 4.6e+18,  # alpha
        "daughter": None
    },
    "76Ge": {
        "name": "Germanium 76",
        "half_life": 1.88e+21,  # double beta-
        "daughter": None
    },
    "174Hf": {
        "name": "Hafnium 174",
        "half_life": 2e+15,  # alpha
        "daughter": None
    },
    "115In": {
        "name": "Indium 115",
        "half_life": 4.41e+14,  # beta-
        "daughter": None
    },
    "78Kr": {
        "name": "Krypton 78",
        "half_life": 9.2e+21,  # double electron capture; NUBASE2020 gives only > 1.1e20 y
        "daughter": None
    },
    "100Mo": {
        "name": "Molybdenum 100",
        "half_life": 7.07e+18,  # double beta-
        "daughter": None
    },
    "144Nd": {
        "name": "Neodymium 144",
        "half_life": 2.29e+15,  # alpha
        "daughter": None
    },
    "150Nd": {
        "name": "Neodymium 150",
        "half_life": 9.3e+18,  # double beta-
        "daughter": None
    },
    "186Os": {
        "name": "Osmium 186",
        "half_life": 2e+15,  # alpha
        "daughter": None
    },
    "148Sm": {
        "name": "Samarium 148",
        "half_life": 6.3e+15,  # alpha
        "daughter": None
    },
    "82Se": {
        "name": "Selenium 82",
        "half_life": 8.76e+19,  # double beta-
        "daughter": None
    },
    "128Te": {
        "name": "Tellurium 128",
        "half_life": 2.25e+24,  # double beta-
        "daughter": None
    },
    "130Te": {
        "name": "Tellurium 130",
        "half_life": 7.91e+20,  # double beta-
        "daughter": None
    },
    "180W": {
        "name": "Tungsten 180",
        "half_life": 1.59e+18,  # alpha
        "daughter": None
    },
    "50V": {
        "name": "Vanadium 50",
        "half_life": 2.71e+17,  # electron capture
        "daughter": None
    },
    "124Xe": {
        "name": "Xenon 124",
        "half_life": 1.8e+22,  # double electron capture (XENON1T, 2019); NUBASE2020 gives only a lower limit
        "daughter": None
    },
    "136Xe": {
        "name": "Xenon 136",
        "half_life": 2.18e+21,  # double beta-
        "daughter": None
    },
    "96Zr": {
        "name": "Zirconium 96",
        "half_life": 2.34e+19,  # double beta-
        "daughter": None
    },
    "238U": {
        "name": "Uranium 238",
        "half_life": 4.463e+09,  # years, alpha
        "daughter": "234Th",
        "daughter_half_life": 0.066  # years (Th-234: 24.1 days)
    },
    "235U": {
        "name": "Uranium 235",
        "half_life": 7.04e+08,  # years, alpha
        "daughter": "231Th",
        "daughter_half_life": 0.002911  # years (Th-231: 25.5 hours)
    },
    "232Th": {
        "name": "Thorium 232",
        "half_life": 1.4e+10,  # years, alpha
        "daughter": "228Ra",
        "daughter_half_life": 5.75  # years (Ra-228: 5.75 years)
    }
}
