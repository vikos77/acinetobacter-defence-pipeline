#!/bin/bash

# Script to build all containers for the containerized pipeline

echo "Building containerized Acinetobacter defense analysis pipeline..."
echo "This will build all required Docker containers."

# Get current user/group IDs for proper permissions
USER_ID=$(id -u)
GROUP_ID=$(id -g)
echo "Building with USER_ID=$USER_ID and GROUP_ID=$GROUP_ID for proper permissions"

# Function to build a container with error checking
build_container() {
    local dockerfile=$1
    local tag=$2
    local context=${3:-.}
    
    echo ""
    echo "===================="
    echo "Building $tag..."
    echo "Dockerfile: $dockerfile"
    echo "===================="
    
    if docker build -f "$dockerfile" \
        --build-arg USER_ID=$USER_ID --build-arg GROUP_ID=$GROUP_ID \
        -t "$tag" "$context"; then
        echo "✅ Successfully built $tag"
    else
        echo "❌ Failed to build $tag"
        return 1
    fi
}

# Function to build a container with target with error checking
build_container_with_target() {
    local dockerfile=$1
    local tag=$2
    local target=$3
    local context=${4:-.}
    
    echo ""
    echo "===================="
    echo "Building $tag..."
    echo "Dockerfile: $dockerfile"
    echo "Target: $target"
    echo "===================="
    
    if docker build -f "$dockerfile" --target "$target" \
        --build-arg USER_ID=$USER_ID --build-arg GROUP_ID=$GROUP_ID \
        -t "$tag" "$context"; then
        echo "✅ Successfully built $tag"
    else
        echo "❌ Failed to build $tag"
        return 1
    fi
}

# Build all containers
echo "Starting container build process..."

# Build in dependency order (some containers don't depend on others, but this is a safe order)
build_container "containers/Dockerfile.edirect" "acinetobacter-edirect:latest" || exit 1
build_container_with_target "containers/Dockerfile.base" "acinetobacter-defensefinder:dev" "development" || exit 1
build_container "containers/Dockerfile.resfinder" "acinetobacter-resfinder:latest" || exit 1
build_container "containers/Dockerfile.padloc" "acinetobacter-padloc:latest" || exit 1
build_container "containers/Dockerfile.crisprcasfinder" "acinetobacter-crisprcasfinder:latest" || exit 1
build_container "containers/Dockerfile.blast" "acinetobacter-blast:latest" || exit 1
build_container "containers/Dockerfile.analysis" "acinetobacter-analysis:latest" || exit 1
build_container "containers/Dockerfile.snakemake" "acinetobacter-snakemake:latest" || exit 1

echo ""
echo "========================================="
echo "✅ All containers built successfully!"
echo "========================================="
echo ""
echo "Built containers:"
docker images | grep "acinetobacter-"
echo ""
echo "Next steps:"
echo "1. Start the services: docker-compose -f docker-compose.containerized.yml up -d"
echo "2. Run the pipeline: docker-compose -f docker-compose.containerized.yml exec snakemake snakemake -f Snakefile.containerized --cores 4"
echo ""
echo "Or use the run_pipeline.sh script for convenience."