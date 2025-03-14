import requests

# Base URL where the API is exposed
BASE_URL = "http://localhost:5051"

def test_health():
    """Test the health check endpoint."""
    url = f"{BASE_URL}/health"
    response = requests.get(url)
    print("Health endpoint:")
    print("Status code:", response.status_code)
    print("Response:", response.json())
    print("-" * 40)

def test_calculate():
    """Test the /calculate endpoint."""
    url = f"{BASE_URL}/calculate"
    payload = {
        "composition": [
            {"fluid": "CO2", "fraction": 0.9},
            {"fluid": "NITROGEN", "fraction": 0.1}
        ],
        "pressure_range": {
            "from": 10,    # in bar
            "to": 300
        },
        "temperature_range": {
            "from": -60,   # in °C
            "to": 150
        },
        "pressure_resolution": 2,
        "temperature_resolution": 5,
        "properties": [
            "enthalpy",
            "vapor_fraction",
            "phase",
            "critical_temperature",
            "critical_pressure",
            "critical_density"
        ],
        "units_system": "SI"
    }
    response = requests.post(url, json=payload)
    print("Calculate endpoint:")
    print("Status code:", response.status_code)
    data = response.json()
    if "results" in data:
        print("Number of results:", len(data["results"]))
        print("First result:", data["results"][0])
    else:
        print("Response:", data)
    print("-" * 40)

def test_ph_flash():
    """Test the ph_flash endpoint.
    
    Adjust the URL below if your blueprint is registered differently.
    For example, if you set a url_prefix '/ph_flash' when creating the blueprint,
    then use:
    
      url = f"{BASE_URL}/ph_flash"
    
    Otherwise, if no prefix is set and the endpoint is registered at '/', change accordingly.
    """

if __name__ == '__main__':
    test_health()
    test_calculate()
