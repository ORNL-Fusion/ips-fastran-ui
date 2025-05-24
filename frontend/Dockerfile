from node:latest

workdir /app

run npm create vite@latest ahoy-app -- --template react

workdir /app/ahoy-app

run npm install

copy react_template/src/App.jsx src/App.jsx 

expose 3000

#cmd ["npm", "run", "dev", "--", "--host"]
cmd ["npm", "run", "dev", "--", "--port", "3000", "--host" ]

