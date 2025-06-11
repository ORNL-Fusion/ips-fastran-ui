from fastapi import FastAPI
from fastapi.responses import Response
from fastapi.middleware.cors import CORSMiddleware
import matplotlib.pyplot as plt
import numpy as np
import io
import base64

app = FastAPI()

origins = [ "http://localhost:3000" ]

app.add_middleware( \
    CORSMiddleware,\
    allow_origins=origins,\
    allow_credentials=True,\
    allow_methods=["*"],\
    allow_headers=["*"],\
    )

@app.get("/", response_class=Response)
async def read_root():

# the below would need to be a request sent to sfapi - if it's too big 
# golden rule: if you would not run this in a login node don't run it in SPIN

  x = np.linspace(0, 10, 100)
  y = np.sin(x)

  plt.figure()
  plt.plot(x, y, label='Sine wave')
  plt.title('Sine wave')
  plt.xlabel('x')
  plt.ylabel('sin(x)')
  plt.legend()

  img_string = io.StringIO()
  plt.savefig(img_string, format='svg')
  plt.close()

  img_content = img_string.getvalue()

  return Response( content=img_content, media_type="image/svg+xml" )
