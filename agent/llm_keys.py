"""Load API keys from environment or project .env at runtime. Keys are never printed, logged, or written."""
import os

def _dotenv():
    vals = {}
    p = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".env")
    if os.path.exists(p):
        for line in open(p, encoding="utf-8"):
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                vals[k.strip()] = v.strip().strip('"').strip("'")
    return vals

def get_key(provider):
    env = {**_dotenv(), **os.environ}
    names = {"gemini": ["GEMINI_API_KEY", "GOOGLE_API_KEY", "GeminiAPI"],
             "deepseek": ["DEEPSEEK_API_KEY"]}[provider]
    for n in names:
        if env.get(n):
            return env[n]
    return None

ENDPOINTS = {"gemini": "https://generativelanguage.googleapis.com/v1beta/openai/",
             "deepseek": "https://api.deepseek.com"}
