# 📂 Structure

## `pyThermoModels/`

- **`__init__.py`**: Initializes the `pyThermoModels` package.
- **`app.py`**: Contains the main entry points for initializing and using the package.
- **`configs/`**: Stores configuration files and constants used throughout the project.
  - `config.py`: Application configuration settings.
  - `constants.py`: Constants such as gas constants and predefined parameters.
  - `setting.py`: Metadata about the package (e.g., version, description).
- **`docs/`**: Contains modules for thermodynamic calculations.
  - `eoscore.py`: Core functionality for equations of state.
  - `eosmanager.py`: Manages EOS-related calculations.
  - `eosmodels.py`: Defines various EOS models.
  - `eosutils.py`: Utility functions for EOS calculations.
  - `thermodb.py`: Handles thermodynamic database operations.
  - `thermomodelcore.py`: Core class for thermodynamic model calculations.
  - `component.py`: Defines components and their thermodynamic properties.
- **`plugin/`**: Contains plugins and reference files.
- **`utils/`**: Utility functions for common operations.

## `test/`

- Contains test scripts and data for validating the functionality of the package.
  - `activity-1.py`, `activity-2.py`, etc.: Test scripts for activity calculations.
  - `thermodb/`: Stores test data for thermodynamic database operations.
