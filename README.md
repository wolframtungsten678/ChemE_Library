# ChemE_Library

## Build Zig

```sh
zig build
```

## Setting up Python

1. Install `python3` developer tools

   ```sh
   sudo apt install python3-dev
   ```

1. Create Python virtual environment

   ```sh
   python3 -m venv .venv
   ```

1. Activate the virtual environment

   ```sh
   source .venv/bin/activate
   ```

1. Load the requirements file

   ```sh
   pip install -r requirements.txt
   ```

## Updating the Requirements.txt

```sh
pip freeze > requirements.txt
```

## Testing

### Adding a file for testing

In `root.zig` add this to get your tests to be detected.

```zig
test {
    ...
    _ = @import("path/to/file.zig");
    ...
}
```

### Minimal Output

```sh
zig build test
```

### Check how many tests runs

```sh
zig build test --summary new
```
