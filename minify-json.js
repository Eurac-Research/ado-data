const fs = require('fs');
const path = require('path');
const { execSync } = require('child_process');

// Function to minify a JSON file
function minifyJsonFile(filePath) {
  try {
    const content = fs.readFileSync(filePath, 'utf8');
    const parsed = JSON.parse(content);
    const minified = JSON.stringify(parsed);

    // Only write if the content is different (to avoid unnecessary file changes)
    if (content !== minified) {
      fs.writeFileSync(filePath, minified);
      return true; // File was modified
    }
    return false; // File was already minified
  } catch (error) {
    console.error(`Error minifying ${filePath}:`, error.message);
    return false;
  }
}

// Function to recursively find all JSON and GeoJSON files
function findJsonFiles(dir) {
  const files = [];

  function traverse(currentDir) {
    const entries = fs.readdirSync(currentDir);

    for (const entry of entries) {
      const fullPath = path.join(currentDir, entry);
      const stat = fs.statSync(fullPath);

      if (stat.isDirectory()) {
        traverse(fullPath);
      } else if (entry.endsWith('.json') || entry.endsWith('.geojson')) {
        files.push(fullPath);
      }
    }
  }

  traverse(dir);
  return files;
}

// Main execution
console.log('Starting JSON/GeoJSON minification...');

const startTime = Date.now();
const allFiles = findJsonFiles('.');
console.log(`Found ${allFiles.length} JSON/GeoJSON files to process`);

let minifiedCount = 0;
let errorCount = 0;

// Process files in batches to avoid memory issues
const batchSize = 50;
for (let i = 0; i < allFiles.length; i += batchSize) {
  const batch = allFiles.slice(i, i + batchSize);

  console.log(`Processing batch ${Math.floor(i / batchSize) + 1}/${Math.ceil(allFiles.length / batchSize)} (${batch.length} files)`);

  for (const file of batch) {
    try {
      const wasMinified = minifyJsonFile(file);
      if (wasMinified) {
        minifiedCount++;
        console.log(`✓ Minified: ${file}`);
      } else {
        console.log(`- Already minified: ${file}`);
      }
    } catch (error) {
      errorCount++;
      console.error(`✗ Error processing ${file}:`, error.message);
    }
  }
}

const endTime = Date.now();
const duration = ((endTime - startTime) / 1000).toFixed(2);

console.log('\n=== Minification Summary ===');
console.log(`Total files processed: ${allFiles.length}`);
console.log(`Files minified: ${minifiedCount}`);
console.log(`Files already minified: ${allFiles.length - minifiedCount - errorCount}`);
console.log(`Errors: ${errorCount}`);
console.log(`Duration: ${duration} seconds`);

// Show some file size statistics
console.log('\n=== File Size Analysis ===');
try {
  const sizeOutput = execSync('du -sh json/', { encoding: 'utf8' });
  console.log(`Total JSON directory size: ${sizeOutput.trim()}`);
} catch (error) {
  console.log('Could not calculate directory size');
}

console.log('\nMinification completed!');
