from setuptools import setup, find_namespace_packages, find_packages

setup(
        name='watertap3',
        url='https://github.com/NREL/WaterTAP3',
        description="WaterTAP3 modeling library",
    keywords="water systems, chemical engineering, process modeling, filtration, desalination, nawi",
        version='0.0.1.dev0',
        packages=find_packages(include=("watertap3*")),
    # install_requires=["pyomo>=6.6.1", "idaes-pse", "numpy", "pandas", "scikit-learn"]
    install_requires=["watertap==0.11", "ipykernel", "scikit-learn"],
        author='WaterTAP3 contributors',
        python_requires=">=3.8",
        author_email='kurban.sitterley@nrel.gov',
        )