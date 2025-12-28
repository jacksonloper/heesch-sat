#!/usr/bin/env node
/**
 * Build script to concatenate all JSON witness files into a single JSONL file
 * for the website to consume.
 */

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const renderingsDir = path.join(__dirname, '../../renderings');
const outputDir = path.join(__dirname, '../public/data');
const outputFile = path.join(outputDir, 'witnesses.jsonl');

// Ensure output directory exists
if (!fs.existsSync(outputDir)) {
  fs.mkdirSync(outputDir, { recursive: true });
}

// Get all JSON files from renderings directory
const jsonFiles = fs.readdirSync(renderingsDir)
  .filter(f => f.endsWith('.json'))
  .sort();

console.log(`Found ${jsonFiles.length} JSON files in ${renderingsDir}`);

// Write JSONL file
const outputStream = fs.createWriteStream(outputFile);

for (const file of jsonFiles) {
  const filePath = path.join(renderingsDir, file);
  const content = fs.readFileSync(filePath, 'utf-8');
  // Parse and re-stringify to ensure single line
  const data = JSON.parse(content);
  outputStream.write(JSON.stringify(data) + '\n');
  console.log(`  Added: ${file}`);
}

outputStream.end();
console.log(`\nWritten to ${outputFile}`);
