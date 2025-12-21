#!/bin/bash
# Complete Podman rootless setup for GAPMananas pipeline
# Run this ONCE on your host system before opening the dev container

set -e

echo "=== Installing Podman and dependencies ==="
sudo apt update
sudo apt install -y \
    podman \
    podman-docker \
    slirp4netns \
    fuse-overlayfs \
    uidmap

echo "=== Configuring subuid/subgid for rootless Podman ==="
# Ensure user has subuid/subgid mappings
if ! grep -q "^${USER}:" /etc/subuid; then
    echo "${USER}:100000:65536" | sudo tee -a /etc/subuid
fi
if ! grep -q "^${USER}:" /etc/subgid; then
    echo "${USER}:100000:65536" | sudo tee -a /etc/subgid
fi

echo "=== Enabling Podman user socket service ==="
# Enable lingering (allows services to run when not logged in)
sudo loginctl enable-linger ${USER}

# Enable and start Podman socket
systemctl --user enable podman.socket
systemctl --user start podman.socket

# Verify socket is running
systemctl --user status podman.socket --no-pager

echo "=== Setting up Podman configuration ==="
mkdir -p ~/.config/containers

# Configure storage to use fuse-overlayfs for better compatibility
cat > ~/.config/containers/storage.conf << EOF
[storage]
driver = "overlay"
runroot = "/run/user/$(id -u)/containers"
graphroot = "${HOME}/.local/share/containers/storage"

[storage.options]
mount_program = "/usr/bin/fuse-overlayfs"
EOF

# Configure registries
cat > ~/.config/containers/registries.conf << 'EOF'
unqualified-search-registries = ["docker.io", "quay.io", "mcr.microsoft.com", "ghcr.io"]
EOF

echo "=== Verifying Podman setup ==="
podman version
echo ""
echo "Socket location:"
ls -la /run/user/$(id -u)/podman/podman.sock

echo ""
echo "=== Testing Podman ==="
podman run --rm hello-world

echo ""
echo "========================================"
echo "✅ Podman rootless setup complete!"
echo "========================================"
echo ""
echo "Podman socket: unix:///run/user/$(id -u)/podman/podman.sock"
echo ""
echo "Next steps:"
echo "1. Open VS Code"
echo "2. Open the GAPMananas folder"
echo "3. Click 'Reopen in Container' when prompted"
echo "4. Wait for post-create.sh to install R/Python packages"
echo ""
