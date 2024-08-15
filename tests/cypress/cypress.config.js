/**
 * Cypress Configuration File
 *
 * This file specifies the setup for running end-to-end tests with Cypress.
 *
 * Base URL is set to 'http://localhost:8902' to match the Flask application port.
 * The spec pattern is set to locate test files in the 'tests/cypress/e2e' directory.
 *
 * To run Cypress with this configuration, use the following command:
 * npx cypress open --config-file tests/cypress/cypress.config.js
 * or
 * npx cypress run --config-file tests/cypress/cypress.config.js
 */


const { defineConfig } = require('cypress')

module.exports = defineConfig({
  e2e: {
    setupNodeEvents(on, config) {},
    baseUrl: 'http://localhost:8902',
    "supportFile": false,
    specPattern: 'e2e/**/*.cy.{js,jsx,ts,tsx}',
  },
})
