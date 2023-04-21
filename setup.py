from setuptools import setup
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="target_spotter",
    version="0.1.0",
    packages=["target_spotter"],
    python_requires=">=3.8",
    package_data={"": ["LICENSE", "*.md", "*.ipynb", "*.yml"]},
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "scikit-learn",
        "joblib",
        "tqdm",
        "statsmodels",
        "glimix-core",
        "numpy-sugar",
        "sqlalchemy",
    ],
    author="Miquel Anglada Girotto",
    author_email="miquel.anglada@crg.eu",
    description="Systematic prioritization of splicing targets to treat cancer.",
    url="https://github.com/MiqG/target_spotter",
    project_urls={"Issues": "https://github.com/MiqG/target_spotter/issues"},
    long_description=long_description,
    long_description_content_type="text/markdown",
)
