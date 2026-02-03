# Docker Build and Release Guide

## Building Docker Images

### Build for Current Version (0.3.4.2)

```bash
# Build the image (run from qtl/ root directory, not qtl/docker/)
docker build -t "saigeqtl" -f docker/Dockerfile .

# Tag with both version number and latest
docker tag saigeqtl wzhou88/saigeqtl:0.3.4.2
docker tag saigeqtl wzhou88/saigeqtl:latest

# Push both tags to Docker Hub
docker push wzhou88/saigeqtl:0.3.4.2
docker push wzhou88/saigeqtl:latest
```

### Build for New Version (Template)

```bash
# Replace X.Y.Z with your version number
VERSION=X.Y.Z

# Build the image
docker build -t "saigeqtl" -f docker/Dockerfile .

# Tag with both version number and latest
docker tag saigeqtl wzhou88/saigeqtl:$VERSION
docker tag saigeqtl wzhou88/saigeqtl:latest

# Push both tags to Docker Hub  
docker push wzhou88/saigeqtl:$VERSION
docker push wzhou88/saigeqtl:latest
```

## Testing Docker Images

### Test the Built Image

```bash
# Test version-specific tag
docker run --rm wzhou88/saigeqtl:0.3.4.2 step1_fitNULLGLMM_qtl.R --help

# Test latest tag
docker run --rm wzhou88/saigeqtl:latest step1_fitNULLGLMM_qtl.R --help

# Test with example data
docker run --rm wzhou88/saigeqtl:latest pixi run R -e "library(SAIGEQTL); packageVersion('SAIGEQTL')"
```

## Available Tags

- `wzhou88/saigeqtl:latest` - Always points to the most recent stable release
- `wzhou88/saigeqtl:0.3.4.2` - Specific version 0.3.4.2 (current stable, adds --solverMethod flag)
- `wzhou88/saigeqtl:0.3.4` - Previous version 0.3.4
- `wzhou88/saigeqtl:0.3.2` - Previous version 0.3.2

## Notes

- **Build Context**: Always run `docker build` from the qtl/ root directory, not from qtl/docker/
- **Dockerfile Location**: Use `-f docker/Dockerfile` to specify the Dockerfile path
- **Multi-platform**: The Dockerfile uses Ubuntu 20.04 base and should work on most platforms
- **Size**: Final image is approximately ~2GB including R, dependencies, and SAIGE-QTL

## Automated Builds

For automated builds via GitHub Actions, the version tagging happens automatically when you run:
```bash
./scripts/update_version.sh X.Y.Z
```

This triggers CI/CD that builds and pushes both `wzhou88/saigeqtl:X.Y.Z` and `wzhou88/saigeqtl:latest`.
