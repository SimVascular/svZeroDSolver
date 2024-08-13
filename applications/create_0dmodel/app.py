from flask import Flask, jsonify, request

app = Flask(__name__)

@app.route('/')
def index():
    return app.send_static_file('index.html')

if __name__ == "__main__":
    # Generate a random port number between 5000 and 9999
    app.run(debug=True, port=8902)
