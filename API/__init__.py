from flask import Flask
from flask_cors import CORS

def create_app():
    """Application factory function."""
    app = Flask(__name__)
    CORS(app)
    
    # Initialize REFPROP
    from API.refprop_setup import RP
    
    # Register all blueprints
    from API.endpoints import blueprints
    for blueprint in blueprints:
        app.register_blueprint(blueprint)
    
    @app.route('/health', methods=['GET'])
    def health_check():
        """Simple health check endpoint."""
        return {'status': 'ok', 'message': 'REFPROP API is running'}
    
    return app