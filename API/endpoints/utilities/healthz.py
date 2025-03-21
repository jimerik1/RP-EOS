from flask import Blueprint, jsonify, current_app, request
import sys
import platform
import datetime
import os
import traceback
from pathlib import Path
from API.refprop_setup import RP
from API.endpoints.utilities.available_fluids import get_fluids_directory

# Create a blueprint for the healthz endpoint
healthz_bp = Blueprint('healthz', __name__)

@healthz_bp.route('/healthz', methods=['GET'])
def healthz():
    """
    Health check endpoint.
    Returns detailed status information about the API.
    
    Query parameters:
    - verbose: If 'true', returns additional system information (default: false)
    """
    try:
        # Parse query parameters
        verbose = request.args.get('verbose', 'false').lower() == 'true'
        
        # Check if REFPROP is initialized by calling a simple function
        wmm = RP.WMOLdll([1.0])
        
        # Get basic system information
        system_info = {
            'timestamp': datetime.datetime.now().isoformat(),
            'api_version': '1.0.0',  # You can update this manually or add version tracking
            'uptime': get_uptime()
        }
        
        # Add more detailed system info if verbose is true
        if verbose:
            system_info.update({
                'python_version': sys.version,
                'platform': platform.platform(),
                'hostname': platform.node(),
                'processor': platform.processor(),
                'memory': get_memory_info()
            })
        
        # Check available endpoints
        endpoints = list_endpoints()
        
        # Check for FLUIDS directory and count fluid files
        fluids_info = get_fluids_info()
        
        # Prepare response
        response = {
            'status': 'ok',
            'message': 'Span-Wagner EOS API is healthy',
            'refprop_status': 'initialized',
            'system_info': system_info,
            'api_info': {
                'endpoints': endpoints,
                'fluids': fluids_info
            }
        }
        
        return jsonify(response)
    except Exception as e:
        traceback.print_exc()
        return jsonify({
            'status': 'error',
            'message': f'Health check failed: {str(e)}',
            'refprop_status': 'not initialized',
            'timestamp': datetime.datetime.now().isoformat()
        }), 500

def get_uptime():
    """Get system uptime if available"""
    try:
        if sys.platform.startswith('linux'):
            with open('/proc/uptime', 'r') as f:
                uptime_seconds = float(f.readline().split()[0])
                return format_uptime(uptime_seconds)
        return "Unavailable"
    except Exception:
        return "Unavailable"

def format_uptime(seconds):
    """Format uptime in a human-readable form"""
    days, remainder = divmod(seconds, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    parts = []
    if days > 0:
        parts.append(f"{int(days)} days")
    if hours > 0:
        parts.append(f"{int(hours)} hours")
    if minutes > 0:
        parts.append(f"{int(minutes)} minutes")
    if seconds > 0 or not parts:
        parts.append(f"{int(seconds)} seconds")
    
    return ", ".join(parts)

def get_memory_info():
    """Get memory information if available"""
    try:
        if sys.platform.startswith('linux'):
            with open('/proc/meminfo', 'r') as f:
                meminfo = dict((i.split()[0].rstrip(':'), int(i.split()[1])) 
                              for i in f.readlines())
                return {
                    'total': f"{meminfo['MemTotal'] // 1024} MB",
                    'free': f"{meminfo['MemFree'] // 1024} MB",
                    'available': f"{meminfo.get('MemAvailable', meminfo['MemFree']) // 1024} MB"
                }
        return "Unavailable"
    except Exception:
        return "Unavailable"

def list_endpoints():
    """List available API endpoints"""
    try:
        app = current_app
        endpoints = []
        
        for rule in app.url_map.iter_rules():
            if not rule.endpoint.startswith('static'):
                endpoints.append({
                    'path': str(rule),
                    'methods': [method for method in rule.methods if method not in ('HEAD', 'OPTIONS')]
                })
        
        return endpoints
    except Exception:
        return "Unavailable"

def get_fluids_info():
    """Get information about available fluids"""
    try:
        fluids_dir = get_fluids_directory()
        
        if fluids_dir.exists():
            fld_files = list(fluids_dir.glob('*.FLD'))
            return {
                'directory': str(fluids_dir),
                'count': len(fld_files),
                'status': 'available',
                'example_fluids': [os.path.splitext(f.name)[0] for f in fld_files[:5]]
            }
        else:
            return {
                'directory': str(fluids_dir),
                'status': 'not found'
            }
    except Exception:
        return {
            'status': 'error'
        }