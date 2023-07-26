from flask import Flask, request, render_template

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('upload.html')

@app.route('/process', methods=['POST'])
def process_file():
    uploaded_file = request.files['upload_file']
    file_content = uploaded_file.read().decode('utf-8')

    # Run your Python script on the file content and get the result
    # Replace the following line with your custom script logic
    processed_result = process_file_content(file_content)

    return processed_result

def process_file_content(content):
    # Custom Python script logic to process the file content
    # In this example, we simply return the content as is
    #print(content)
    return content

if __name__ == '__main__':
    app.run()

