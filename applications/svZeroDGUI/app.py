"""
app.py

This script sets up a basic Flask web server for serving static files and handling
HTTP requests. It includes:

- Initialization of the Flask application.
- A route for serving the main HTML file (`index.html`) when accessing the root URL ('/').
- Configuration to run the application with debugging enabled on port 8902.
"""

from flask import Flask, jsonify, request

app = Flask(__name__)

@app.route('/')
def index():
    return app.send_static_file('index.html')

if __name__ == "__main__":
    # Generate a random port number between 5000 and 9999
    app.run(debug=True, port=8902)
