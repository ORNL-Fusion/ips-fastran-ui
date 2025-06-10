from python:slim-bookworm

copy requirements.txt .

run pip install --no-cache-dir -r requirements.txt

copy app.py .

expose 8000

cmd ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000", "--reload"]
