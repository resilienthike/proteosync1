import sys

print("Python:", sys.executable)

# Example: import from your package (reflects edits in src/vsx)
try:
    import vsx

    print("vsx version:", getattr(vsx, "__version__", "unknown"))
except Exception as e:
    print("could not import vsx:", e)
