from flask import Flask, request

app = Flask(__name__)

@app.route('/process', methods=['POST'])
def process_file():
    uploaded_file = request.files['upload_file']
    # Perform the processing on the uploaded file (e.g., read content, apply logic)
    result = process_uploaded_file(uploaded_file)
    return result

def process_uploaded_file(file):
    # Perform your custom processing logic here
    # Example: Read file content and return it as a result
    content = file.read().decode('utf-8')
    return content

if __name__ == '__main__':
    app.run()

