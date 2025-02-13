from pathlib import Path
import yaml
from typing import Dict

class DefaultsLoader:
    """Loader for YAML default configurations."""
    
    def __init__(self):
        self._defaults_dir = Path(__file__).parent
        self._cache: Dict[str, dict] = {}
    
    def load(self, name: str) -> dict:
        # Ensure .yaml extension
        if not name.endswith('.yaml'):
            name = f"{name}.yaml"
            
        # Return cached config if available
        if name in self._cache:
            return self._cache[name]
            
        # Load and cache new config
        file_path = self._defaults_dir / name
        if not file_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {file_path}")
            
        with open(file_path) as f:
            config = yaml.safe_load(f)
            self._cache[name] = config
            return config

# Create a singleton instance
loader = DefaultsLoader()