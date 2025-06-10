docker run -it --rm \
-v $(pwd):/host \
--entrypoint sh \
--user $(id -u):$(id -g) \
react_example
#docker run -d -p 3000:3000 --name react-app --link fastapi:fastapi react_example
#npm init vite@latest frontend -- --template react
