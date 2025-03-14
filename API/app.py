from API import create_app

# Create the Flask application
app = create_app()

if __name__ == '__main__':
    # Bind to 0.0.0.0 so Docker can map the port externally
    app.run(host='0.0.0.0', debug=False, use_reloader=False)