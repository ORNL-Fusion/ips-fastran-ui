#docker build -t fastapi_example .

# new docker build command below, original old above

# registry/project/namespace/image_name:tag <--- namespace is different than spin namespace
# registry.nersc.gov has a web interface; the namespace can be modified there

docker build -t registry.nersc.gov/atom/ips_fastran_test/fastapi_example:latest .
