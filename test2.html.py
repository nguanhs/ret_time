from flask import Flask, request

app = Flask(__name__)

@app.route('/process', methods=['POST'])
def process_file():
    uploaded_file = request.files['upload_file']
    content = uploaded_file.read().decode('utf-8')
    return content

if __name__ == '__main__':
    app.run()

