FROM julia:1.8.0

WORKDIR /app

ADD example/ .

RUN julia -e 'using Pkg;Pkg.add(PackageSpec(url="https://github.com/KrainskiL/CGE.jl"));'