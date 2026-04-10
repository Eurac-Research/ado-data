const INVALID_JSON_VALUE_PATTERN = /^[-+]?(?:Infinity|NaN|infitiy)/i;

function findPreviousSignificantChar(content, startIndex) {
  for (let index = startIndex - 1; index >= 0; index -= 1) {
    if (!/\s/.test(content[index])) {
      return content[index];
    }
  }

  return null;
}

function findNextSignificantChar(content, startIndex) {
  for (let index = startIndex; index < content.length; index += 1) {
    if (!/\s/.test(content[index])) {
      return content[index];
    }
  }

  return null;
}

function sanitizeInvalidJsonValues(content) {
  let sanitized = '';
  let replacementCount = 0;
  let inString = false;
  let isEscaped = false;

  for (let index = 0; index < content.length; ) {
    const char = content[index];

    if (inString) {
      sanitized += char;

      if (isEscaped) {
        isEscaped = false;
      } else if (char === '\\') {
        isEscaped = true;
      } else if (char === '"') {
        inString = false;
      }

      index += 1;
      continue;
    }

    if (char === '"') {
      inString = true;
      sanitized += char;
      index += 1;
      continue;
    }

    const match = content.slice(index).match(INVALID_JSON_VALUE_PATTERN);

    if (match) {
      const previousChar = findPreviousSignificantChar(content, index);
      const nextChar = findNextSignificantChar(content, index + match[0].length);
      const isValuePosition = previousChar === ':' || previousChar === ',' || previousChar === '[';
      const hasValidTerminator = nextChar === ',' || nextChar === '}' || nextChar === ']' || nextChar === null;

      if (isValuePosition && hasValidTerminator) {
        sanitized += 'null';
        replacementCount += 1;
        index += match[0].length;
        continue;
      }
    }

    sanitized += char;
    index += 1;
  }

  return { sanitized, replacementCount };
}

function parseJsonWithFallback(content, filePath, logger = console) {
  try {
    return {
      data: JSON.parse(content),
      sanitized: false,
      replacementCount: 0
    };
  } catch (parseError) {
    const { sanitized, replacementCount } = sanitizeInvalidJsonValues(content);

    if (replacementCount === 0) {
      parseError.message = `Invalid JSON in ${filePath}: ${parseError.message}`;
      throw parseError;
    }

    try {
      const data = JSON.parse(sanitized);

      if (logger && typeof logger.warn === 'function') {
        logger.warn(
          `⚠️  Replaced ${replacementCount} invalid JSON value${replacementCount === 1 ? '' : 's'} with null in ${filePath}`
        );
      }

      return {
        data,
        sanitized: true,
        replacementCount
      };
    } catch (sanitizedParseError) {
      sanitizedParseError.message =
        `Invalid JSON in ${filePath} after sanitizing ${replacementCount} invalid value` +
        `${replacementCount === 1 ? '' : 's'}: ${sanitizedParseError.message}`;
      throw sanitizedParseError;
    }
  }
}

module.exports = {
  parseJsonWithFallback,
  sanitizeInvalidJsonValues
};
