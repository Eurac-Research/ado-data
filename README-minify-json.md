# JSON/GeoJSON Minification Script

A secure and efficient Node.js script for minifying JSON and GeoJSON files in the Alpine Drought Observatory (ADO) data repository.

## Features

- ✅ **Secure Processing**: Only processes `.json` and `.geojson` files with built-in security validations
- 🚀 **Efficient**: Batch processing to handle large numbers of files
- 📊 **Smart**: Skips already minified files to avoid unnecessary processing
- 📈 **Detailed Reporting**: Comprehensive progress tracking and statistics
- 🛡️ **Safe**: Input validation, file size limits, and error handling
- 🔍 **Flexible**: Supports dry-run mode and verbose output

## Usage

```bash
# Basic usage - minify all JSON/GeoJSON files
node minify-json.js

# Preview what would be minified (dry run)
node minify-json.js --dry-run

# Show detailed output for each file
node minify-json.js --verbose

# Show help
node minify-json.js --help
```

## Options

| Option | Short | Description |
|--------|-------|-------------|
| `--help` | `-h` | Show help message |
| `--dry-run` | `-d` | Preview changes without modifying files |
| `--verbose` | `-v` | Show detailed output for each file processed |

## Security Features

- **File Type Validation**: Only processes files with `.json` and `.geojson` extensions
- **File Size Limits**: Prevents processing of files larger than 100MB
- **Path Validation**: Ensures processing stays within the intended directory
- **JSON Validation**: Validates JSON structure before writing
- **Safe Directory Traversal**: Skips hidden directories and common build folders
- **Error Handling**: Graceful handling of file system errors

## Output

The script provides detailed progress reporting:

```
🚀 Starting JSON/GeoJSON minification...
📁 Found 688 JSON/GeoJSON files to process

📦 Processing batch 1/14 (50 files)
✅ Minified: json/hydro/CDI-latest.geojson
⏭️  Already minified: json/nuts/CDI-latest-features.json

==================================================
📊 MINIFICATION SUMMARY
==================================================
📁 Total files processed: 688
✅ Files minified: 633
⏭️  Files already minified: 55
❌ Errors: 0
⏱️  Duration: 19.26 seconds
💾 Size reduction: 245.67MB (23.4% smaller)
```

## Configuration

The script includes configurable constants:

```javascript
const CONFIG = {
  SUPPORTED_EXTENSIONS: ['.json', '.geojson'],
  MAX_FILE_SIZE: 100 * 1024 * 1024, // 100MB limit
  BATCH_SIZE: 50,
  ENCODING: 'utf8'
};
```

## Requirements

- Node.js (version 12 or higher)
- Read/write permissions in the target directory

## How It Works

1. **Discovery**: Recursively finds all JSON and GeoJSON files
2. **Validation**: Checks file type, size, and accessibility
3. **Processing**: Parses JSON and removes unnecessary whitespace
4. **Optimization**: Only writes files that have been changed
5. **Reporting**: Provides detailed statistics and progress updates

## Error Handling

The script includes comprehensive error handling:

- Invalid JSON files are skipped with detailed error messages
- File system errors are logged but don't stop the process
- The script exits with appropriate error codes
- Uncaught exceptions are handled gracefully

## Performance

- **Batch Processing**: Processes files in batches of 50 to manage memory usage
- **Skip Unchanged**: Only processes files that need minification
- **Memory Efficient**: Processes files one at a time within batches
- **File Size Limits**: Prevents processing of extremely large files

## License

MIT License - See the main repository LICENSE file for details.

## Contributing

This script is part of the Alpine Drought Observatory project. For contributions or issues, please refer to the main repository guidelines.
