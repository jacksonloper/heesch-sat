#!/usr/bin/env node
/**
 * Build script to fetch all witness data from Modal API and create a static JSONL file
 * for the website to consume at runtime.
 */

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const outputDir = path.join(__dirname, '../public/data');
const outputFile = path.join(outputDir, 'witnesses.jsonl');

const MODAL_API_URL = 'https://hloper--heesch-renderings-web.modal.run/list_full';

async function fetchFromModal() {
  console.log(`Fetching data from Modal API: ${MODAL_API_URL}`);

  const response = await fetch(MODAL_API_URL);
  if (!response.ok) {
    throw new Error(`Failed to fetch from Modal: ${response.status} ${response.statusText}`);
  }

  const data = await response.json();
  return data.polyforms || [];
}

async function main() {
  // Ensure output directory exists
  if (!fs.existsSync(outputDir)) {
    fs.mkdirSync(outputDir, { recursive: true });
  }

  try {
    const polyforms = await fetchFromModal();
    console.log(`Received ${polyforms.length} polyforms from Modal API`);

    // Write JSONL file
    const outputStream = fs.createWriteStream(outputFile);

    for (const polyform of polyforms) {
      outputStream.write(JSON.stringify(polyform) + '\n');
    }

    outputStream.end();
    console.log(`\nWritten to ${outputFile}`);
  } catch (error) {
    console.error('Error fetching from Modal:', error.message);

    // Fallback: try to read from local renderings directory
    const renderingsDir = path.join(__dirname, '../../renderings');
    if (fs.existsSync(renderingsDir)) {
      console.log('Falling back to local renderings directory...');

      const jsonFiles = fs.readdirSync(renderingsDir)
        .filter(f => f.endsWith('.json'))
        .sort();

      console.log(`Found ${jsonFiles.length} JSON files in ${renderingsDir}`);

      const outputStream = fs.createWriteStream(outputFile);

      for (const file of jsonFiles) {
        const filePath = path.join(renderingsDir, file);
        const content = fs.readFileSync(filePath, 'utf-8');
        const data = JSON.parse(content);
        outputStream.write(JSON.stringify(data) + '\n');
        console.log(`  Added: ${file}`);
      }

      outputStream.end();
      console.log(`\nWritten to ${outputFile}`);
    } else {
      console.error('No local renderings directory found. Build failed.');
      process.exit(1);
    }
  }
}

main();
