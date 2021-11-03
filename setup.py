from setuptools import setup

with open('README.md') as f:
    long_description = f.read()

setup(
    name="target_spotter",
    version="0.1.0",
    packages=["target_spotter"],
    python_requires=">=3.8",
    package_data={"": ["LICENSE", "*.md","*.ipynb","*.yml"]},
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "scikit-learn",
        "joblib",
        "tqdm",
        "statsmodels"
    ],
    author="Miquel Anglada Girotto",
    author_email="miquel.anglada@crg.eu",
    description="Systematic prioritization of splicing targets to treat cancer.",
    url="https://github.com/CRG-CNAG/target_spotter",
    project_urls={"Issues": "https://github.com/CRG-CNAG/target_spotter/issues"},
    long_description=long_description,
    long_description_content_type='text/markdown'
)
