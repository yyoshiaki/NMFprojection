import setuptools

import NMFprojection

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="NMFproj", # Replace with your own username
    version=attr: NMFprojection._version.__version__,
    install_requires=[
        "requests",
        "pandas",
        "numpy",
        "scikit-learn",
        "scanpy"
    ],
    entry_points={
        'console_scripts': [
            'NMFproj=NMFprojection:main',
        ],
    },
    author="Yoshiaki Yasumizu",
    author_email="yyasumizu@ifrec.osaka-u.ac.jp",
    description="NMF projection",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yyoshiaki/NMFprojection",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: Attribution-NonCommercial 4.0 International",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)