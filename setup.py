from setuptools import setup, find_packages

# Read the README file for a detailed project description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="chem_inf_widgets",  # The package name
    version="0.1.0",  # Package version
    author="Crtomir Podlipnik",  # Replace with your name
    author_email="crtomir.podlipnik@fkkt.uni-lj.si",  # Replace with your email
    description="Custom chemoinformatics widgets for Orange3",  # Short description
    long_description=long_description,  # Full description from README.md
    long_description_content_type="text/markdown",  # Format of the README file
    url="https://github.com/crtomirp/chem-inf-widgets",  # Project repository URL
    packages=find_packages(),  # Automatically find all packages in the directory
    include_package_data=True,  # Include additional files (e.g., icons)
    install_requires=[
        "setuptools",
        "Orange3>=3.32.0",  # Dependency for Orange3
        "rdkit",  # Example dependency for chemoinformatics
        "PyQtWebEngine", #
        "requests-cache",
        "pdfkit",
        "pubchempy",
        "mordred",
        "xhtml2pdf",
        "torch==2.3.0",
        # Add additional dependencies here
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",  # Project maturity level
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",  # License type
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",  # Supported Python versions
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.8",  # Minimum Python version required
    license="MIT",  # License type
    keywords="chemoinformatics orange3 widgets sdf",  # Relevant keywords
    entry_points={
        "orange.widgets": [
            "Chemoinformatics = chem_inf_widgets",  # Register widgets in Orange3
        ],
    },
    project_urls={
        "Documentation": "https://github.com/crtomirp/chem-inf-widgets#readme",
        "Source": "https://github.com/crtomirp/chem-inf-widgets",
        "Tracker": "https://github.com/crtomirp/chem-inf-widgets/issues",
    },
)
