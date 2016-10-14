try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Code to analyse overexpression studies',
    'author': 'Vito Zanotelli',
    'url': 'bodenmillerlab.org',
    'download_url': 'Where to download it.',
    'author_email': 'vito.zanotelli@gmail.com',
    'version': '0.1',
    'install_requires': ['nose', 'argparse', 'numpy', 
'pandas', 'seaborn', 'scipy', 'configparser2', 'fcm', 'matplotlib'],
    'packages': ['adnet'],
    'scripts': ['adnet/adnet_analysis.py'],
    'name': 'ADNET'
}

setup(**config)
