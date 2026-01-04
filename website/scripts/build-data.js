#!/usr/bin/env node
/* global process */
/**
 * Build script to read all witness data from local renderings directory
 * and create a static JSONL file for the website to consume at runtime.
 */

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const outputDir = path.join(__dirname, '../public/data');
const outputFile = path.join(outputDir, 'witnesses.jsonl');
const renderingsDir = path.join(__dirname, '../../renderings');

function main() {
  // Ensure output directory exists
  if (!fs.existsSync(outputDir)) {
    fs.mkdirSync(outputDir, { recursive: true });
  }

  if (!fs.existsSync(renderingsDir)) {
    console.error(`Error: renderings directory not found at ${renderingsDir}`);
    process.exit(1);
  }

  // Read all JSON files from local renderings directory
  // Exclude search_* files which are summary/metadata files, not polyform data
  const jsonFiles = fs.readdirSync(renderingsDir)
    .filter(f => f.endsWith('.json') && !f.startsWith('search_'))
    .sort();

  console.log(`Found ${jsonFiles.length} JSON files in ${renderingsDir}`);

  const outputStream = fs.createWriteStream(outputFile);
  let count = 0;

  for (const file of jsonFiles) {
    const filePath = path.join(renderingsDir, file);
    const content = fs.readFileSync(filePath, 'utf-8');
    const data = JSON.parse(content);
    
    // Only include files that have required polyform fields
    if (data.coordinates && data.grid_type) {
      outputStream.write(JSON.stringify(data) + '\n');
      count++;
      console.log(`  Added: ${file}`);
    } else {
      console.log(`  Skipped (not a polyform): ${file}`);
    }
  }

  outputStream.end();
  console.log(`\nWritten ${count} polyforms to ${outputFile}`);
}

main();
