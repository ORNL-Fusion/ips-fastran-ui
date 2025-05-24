# ips-fastran-ui

Description:

  This repo contains two containers, one in directory frontend, one in directory backend.

  The backend container runs a FastAPI program that generates a sine wave with matplotlib,
  saves it as svg xml, and returns it at localhost:8000

  The frontend container runs a react.js program that gets the svg data from localhost:8000
  from the backend and displays it at localhost:3000

  Design justification:

  Why are we plotting data in the backend and sending the svg to the frontend instead of
  sending the data to the frontend and plotting it there? Because this webapp will wrap an
  existing command line workflow in which all the graphs are produced with matplotlib in
  python. We don't want to rewrite all that functionality at least initially.


Instructions:


  docker container ps -a

  (if you see fastapi_container or react_container, run the next 2 steps: )

  docker container kill fastapi_container react_container

  docker container prune

  (then build the containers and run them)

  cd backend;
  bash build.sh;
  bash run.sh;

  cd frontend;
  bash build.sh;
  bash run.sh;



  To see the backend, go to localhost:8000

  To see the frontend, go to localhost:3000

  You should see the same graph

![image](https://github.com/user-attachments/assets/4e3eefed-114b-417f-a654-ca1154df0829)

![image](https://github.com/user-attachments/assets/7f9736b9-e30c-40a5-9f52-117578a84c79)

