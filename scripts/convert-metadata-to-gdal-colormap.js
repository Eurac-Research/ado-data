#!/usr/bin/env node

/**
 * Convert ADO metadata JSON colormap to GDAL color-relief format
 * 
 * Usage: node convert-metadata-to-gdal-colormap.js path/to/VHI.json
 * Output: VHI-colormap.txt (GDAL color-relief format)
 */

const fs = require('fs');
const path = require('path');

function parseRGBA(rgbaString) {
  const match = rgbaString.match(/rgba?\((\d+),\s*(\d+),\s*(\d+)(?:,\s*([\d.]+))?\)/);
  if (!match) return [0, 0, 0, 255];
  return [
    parseInt(match[1]),
    parseInt(match[2]),
    parseInt(match[3]),
    match[4] ? Math.round(parseFloat(match[4]) * 255) : 255
  ];
}

function convertMetadataToGDALColormap(metadataPath, outputDir = null) {
  // Read the metadata JSON
  const metadata = JSON.parse(fs.readFileSync(metadataPath, 'utf8'));
  
  const indexName = metadata.short_name;
  const stops = metadata.colormap?.paint?.['fill-color']?.stops;
  
  if (!stops || !Array.isArray(stops)) {
    throw new Error('Invalid colormap format in metadata');
  }
  
  console.log(`Converting ${indexName} colormap...`);
  console.log(`Found ${stops.length} color stops`);
  
  // Convert to GDAL format
  // Format: value R G B [A]
  const gdalLines = stops.map(([value, color]) => {
    const [r, g, b, a] = parseRGBA(color);
    return `${value} ${r} ${g} ${b} ${a}`;
  });
  
  // Add header comment
  const header = [
    `# GDAL Color Relief for ${indexName}`,
    `# Generated from: ${path.basename(metadataPath)}`,
    `# Date: ${new Date().toISOString()}`,
    `# Format: value R G B A`,
    ''
  ];
  
  const outputContent = [...header, ...gdalLines].join('\n');
  
  // Write output file
  const outputPath = outputDir 
    ? path.join(outputDir, `${indexName}-colormap.txt`)
    : path.join(path.dirname(metadataPath), `${indexName}-colormap.txt`);
  
  fs.writeFileSync(outputPath, outputContent);
  console.log(`✓ Created: ${outputPath}`);
  
  return outputPath;
}

// Main
if (require.main === module) {
  const metadataPath = process.argv[2];
  const outputDir = process.argv[3];
  
  if (!metadataPath) {
    console.error('Usage: node convert-metadata-to-gdal-colormap.js <metadata.json> [outputDir]');
    console.error('Example: node convert-metadata-to-gdal-colormap.js json/nuts/metadata/VHI.json tiffs/colormaps');
    process.exit(1);
  }
  
  if (!fs.existsSync(metadataPath)) {
    console.error(`Error: File not found: ${metadataPath}`);
    process.exit(1);
  }
  
  try {
    convertMetadataToGDALColormap(metadataPath, outputDir);
  } catch (error) {
    console.error('Error:', error.message);
    process.exit(1);
  }
}

module.exports = { convertMetadataToGDALColormap };
