#!/usr/bin/env node

const fs = require('fs');
const path = require('path');
const { parseJsonWithFallback } = require('./json-utils');

function main() {
  const nutsDir = path.join('json', 'nuts');

  if (!fs.existsSync(nutsDir)) {
    console.log('json/nuts directory does not exist');
    process.exit(0);
  }

  const files = fs.readdirSync(nutsDir).filter((file) => file.endsWith('-latest.geojson'));

  if (files.length === 0) {
    console.log('No *-latest.geojson files found in json/nuts directory');
    process.exit(0);
  }

  console.log(`Found ${files.length} *-latest.geojson files to process`);

  let errorCount = 0;

  files.forEach((file) => {
    const filePath = path.join(nutsDir, file);
    console.log(`Processing: ${file}`);

    try {
      const content = fs.readFileSync(filePath, 'utf8');
      const { data: geojsonData } = parseJsonWithFallback(content, filePath, console);

      if (!geojsonData.features || !Array.isArray(geojsonData.features)) {
        console.log(`No features array found in ${file}`);
        return;
      }

      const propertiesOnly = geojsonData.features
        .map((feature) => feature.properties || null)
        .filter((properties) => properties !== null);

      const outputData = {
        features: propertiesOnly
      };

      if (geojsonData.metadata) {
        outputData.metadata = geojsonData.metadata;
      }

      const baseName = path.basename(file, '.geojson');
      const outputFile = path.join(nutsDir, `${baseName}-features.json`);

      fs.writeFileSync(outputFile, JSON.stringify(outputData));

      console.log(
        `Created: ${outputFile} with ${propertiesOnly.length} feature properties` +
          `${geojsonData.metadata ? ' and metadata' : ''} (minified, no geometry)`
      );

    } catch (error) {
      errorCount += 1;
      console.error(`Error processing ${file}: ${error.message}`);
    }
  });

  console.log('Feature extraction completed');
  process.exit(errorCount > 0 ? 1 : 0);
}

main();
