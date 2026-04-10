#!/usr/bin/env node

/**
 * JSON/GeoJSON Minification Script
 * 
 * This script recursively finds and minifies all JSON and GeoJSON files in the current directory
 * and its subdirectories. It removes unnecessary whitespace and formatting while preserving
 * the data structure and content.
 * 
 * Features:
 * - Recursively processes all .json and .geojson files
 * - Batch processing to handle large numbers of files efficiently
 * - Skips already minified files to avoid unnecessary processing
 * - Provides detailed progress reporting and statistics
 * - Safe error handling with detailed error reporting
 * - Input validation and security checks
 * 
 * Usage:
 *   node minify-json.js [options]
 * 
 * Options:
 *   --help, -h     Show this help message
 *   --dry-run, -d  Show what would be minified without making changes
 *   --verbose, -v  Show detailed output for each file processed
 * 
 * Security:
 * - Only processes files with .json and .geojson extensions
 * - Validates JSON structure before writing
 * - Uses safe file system operations
 * - Limits file size processing to prevent memory issues
 * 
 * @author Alpine Drought Observatory Team
 * @version 1.0.0
 * @license MIT
 */

const fs = require('fs');
const path = require('path');
const { execSync } = require('child_process');
const { parseJsonWithFallback } = require('./scripts/json-utils');

// Configuration constants
const CONFIG = {
  SUPPORTED_EXTENSIONS: ['.json', '.geojson'],
  MAX_FILE_SIZE: 100 * 1024 * 1024, // 100MB limit
  BATCH_SIZE: 50,
  ENCODING: 'utf8'
};

// Command line argument parsing
const args = process.argv.slice(2);
const options = {
  help: args.includes('--help') || args.includes('-h'),
  dryRun: args.includes('--dry-run') || args.includes('-d'),
  verbose: args.includes('--verbose') || args.includes('-v')
};

/**
 * Display help information
 */
function showHelp() {
  console.log(`
JSON/GeoJSON Minification Script

Usage: node minify-json.js [options]

Options:
  --help, -h     Show this help message
  --dry-run, -d  Show what would be minified without making changes
  --verbose, -v  Show detailed output for each file processed

Examples:
  node minify-json.js                    # Minify all JSON files
  node minify-json.js --dry-run          # Preview what would be minified
  node minify-json.js --verbose          # Show detailed processing info

Security Features:
- Only processes .json and .geojson files
- File size limit: ${CONFIG.MAX_FILE_SIZE / (1024 * 1024)}MB
- Validates JSON structure before writing
- Safe error handling and reporting
`);
}

/**
 * Validate file for processing
 * @param {string} filePath - Path to the file to validate
 * @returns {boolean} - True if file is safe to process
 */
function validateFile(filePath) {
  try {
    // Check if file exists and is accessible
    const stats = fs.statSync(filePath);
    
    // Security check: Only process regular files
    if (!stats.isFile()) {
      console.warn(`⚠️  Skipping non-file: ${filePath}`);
      return false;
    }
    
    // Security check: File size limit
    if (stats.size > CONFIG.MAX_FILE_SIZE) {
      console.warn(`⚠️  Skipping large file (${(stats.size / (1024 * 1024)).toFixed(2)}MB): ${filePath}`);
      return false;
    }
    
    // Security check: Only process supported extensions
    const ext = path.extname(filePath).toLowerCase();
    if (!CONFIG.SUPPORTED_EXTENSIONS.includes(ext)) {
      console.warn(`⚠️  Skipping unsupported file type: ${filePath}`);
      return false;
    }
    
    return true;
  } catch (error) {
    console.error(`❌ Error validating file ${filePath}:`, error.message);
    return false;
  }
}

/**
 * Minify a JSON file by removing unnecessary whitespace
 * @param {string} filePath - Path to the JSON file to minify
 * @returns {boolean} - True if file was modified, false if already minified or error
 */
function minifyJsonFile(filePath) {
  try {
    // Validate file before processing
    if (!validateFile(filePath)) {
      return false;
    }
    
    // Read file content
    const content = fs.readFileSync(filePath, CONFIG.ENCODING);
    
    // Validate JSON structure
    let parsed;
    try {
      ({ data: parsed } = parseJsonWithFallback(content, filePath, console));
    } catch (parseError) {
      console.error(`❌ ${parseError.message}`);
      return false;
    }
    
    // Minify the JSON
    const minified = JSON.stringify(parsed);
    
    // Only write if the content is different (to avoid unnecessary file changes)
    if (content !== minified) {
      if (options.dryRun) {
        console.log(`🔍 Would minify: ${filePath} (${content.length} → ${minified.length} chars)`);
        return true;
      }
      
      // Write minified content
      fs.writeFileSync(filePath, minified, CONFIG.ENCODING);
      
      if (options.verbose) {
        const savings = ((content.length - minified.length) / content.length * 100).toFixed(1);
        console.log(`✅ Minified: ${filePath} (${content.length} → ${minified.length} chars, ${savings}% smaller)`);
      }
      
      return true; // File was modified
    }
    
    if (options.verbose) {
      console.log(`⏭️  Already minified: ${filePath}`);
    }
    
    return false; // File was already minified
  } catch (error) {
    console.error(`❌ Error minifying ${filePath}:`, error.message);
    return false;
  }
}

/**
 * Recursively find all JSON and GeoJSON files in a directory
 * @param {string} dir - Directory to search
 * @returns {string[]} - Array of file paths
 */
function findJsonFiles(dir) {
  const files = [];
  
  /**
   * Recursively traverse directories
   * @param {string} currentDir - Current directory being processed
   */
  function traverse(currentDir) {
    try {
      // Security check: Ensure we're working within safe bounds
      const resolvedPath = path.resolve(currentDir);
      const basePath = path.resolve(dir);
      
      if (!resolvedPath.startsWith(basePath)) {
        console.warn(`⚠️  Skipping directory outside of base path: ${currentDir}`);
        return;
      }
      
      const entries = fs.readdirSync(currentDir);
      
      for (const entry of entries) {
        const fullPath = path.join(currentDir, entry);
        
        try {
          const stat = fs.statSync(fullPath);
          
          if (stat.isDirectory()) {
            // Skip hidden directories and common non-data directories
            if (!entry.startsWith('.') && !['node_modules', 'dist', 'build'].includes(entry)) {
              traverse(fullPath);
            }
          } else if (stat.isFile()) {
            const ext = path.extname(entry).toLowerCase();
            if (CONFIG.SUPPORTED_EXTENSIONS.includes(ext)) {
              files.push(fullPath);
            }
          }
        } catch (statError) {
          console.warn(`⚠️  Cannot access ${fullPath}:`, statError.message);
        }
      }
    } catch (error) {
      console.error(`❌ Error reading directory ${currentDir}:`, error.message);
    }
  }
  
  traverse(dir);
  return files;
}

/**
 * Main execution function
 */
function main() {
  // Show help if requested
  if (options.help) {
    showHelp();
    process.exit(0);
  }
  
  console.log('🚀 Starting JSON/GeoJSON minification...');
  
  if (options.dryRun) {
    console.log('🔍 DRY RUN MODE: No files will be modified');
  }
  
  const startTime = Date.now();
  const allFiles = findJsonFiles('.');
  
  if (allFiles.length === 0) {
    console.log('ℹ️  No JSON or GeoJSON files found in the current directory');
    process.exit(0);
  }
  
  console.log(`📁 Found ${allFiles.length} JSON/GeoJSON files to process`);
  
  let minifiedCount = 0;
  let errorCount = 0;
  let totalOriginalSize = 0;
  let totalMinifiedSize = 0;
  
  // Process files in batches to avoid memory issues
  const batchSize = CONFIG.BATCH_SIZE;
  const totalBatches = Math.ceil(allFiles.length / batchSize);
  
  for (let i = 0; i < allFiles.length; i += batchSize) {
    const batch = allFiles.slice(i, i + batchSize);
    const batchNumber = Math.floor(i / batchSize) + 1;
    
    console.log(`\n📦 Processing batch ${batchNumber}/${totalBatches} (${batch.length} files)`);
    
    for (const file of batch) {
      try {
        // Track file sizes for statistics
        const originalSize = fs.statSync(file).size;
        totalOriginalSize += originalSize;
        
        const wasMinified = minifyJsonFile(file);
        
        if (wasMinified) {
          minifiedCount++;
          if (!options.dryRun) {
            const minifiedSize = fs.statSync(file).size;
            totalMinifiedSize += minifiedSize;
            
            if (!options.verbose) {
              console.log(`✅ Minified: ${file}`);
            }
          }
        } else {
          totalMinifiedSize += originalSize;
          if (!options.verbose) {
            console.log(`⏭️  Already minified: ${file}`);
          }
        }
      } catch (error) {
        errorCount++;
        console.error(`❌ Error processing ${file}:`, error.message);
      }
    }
  }
  
  const endTime = Date.now();
  const duration = ((endTime - startTime) / 1000).toFixed(2);
  
  // Display summary
  console.log('\n' + '='.repeat(50));
  console.log('📊 MINIFICATION SUMMARY');
  console.log('='.repeat(50));
  console.log(`📁 Total files processed: ${allFiles.length}`);
  console.log(`✅ Files minified: ${minifiedCount}`);
  console.log(`⏭️  Files already minified: ${allFiles.length - minifiedCount - errorCount}`);
  console.log(`❌ Errors: ${errorCount}`);
  console.log(`⏱️  Duration: ${duration} seconds`);
  
  if (!options.dryRun && minifiedCount > 0) {
    const sizeSavings = totalOriginalSize - totalMinifiedSize;
    const percentSavings = ((sizeSavings / totalOriginalSize) * 100).toFixed(1);
    console.log(`💾 Size reduction: ${(sizeSavings / (1024 * 1024)).toFixed(2)}MB (${percentSavings}% smaller)`);
  }
  
  // Show file size statistics
  console.log('\n📈 FILE SIZE ANALYSIS');
  console.log('='.repeat(50));
  try {
    const sizeOutput = execSync('du -sh json/', { encoding: CONFIG.ENCODING });
    console.log(`📂 Total JSON directory size: ${sizeOutput.trim()}`);
  } catch (error) {
    console.log('ℹ️  Could not calculate directory size');
  }
  
  console.log('\n🎉 Minification completed!');
  
  // Exit with appropriate code
  process.exit(errorCount > 0 ? 1 : 0);
}

// Handle uncaught exceptions gracefully
process.on('uncaughtException', (error) => {
  console.error('❌ Uncaught exception:', error.message);
  process.exit(1);
});

process.on('unhandledRejection', (reason, promise) => {
  console.error('❌ Unhandled rejection at:', promise, 'reason:', reason);
  process.exit(1);
});

// Run the main function
if (require.main === module) {
  main();
}
