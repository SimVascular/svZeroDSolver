/**
 * Cypress Configuration File
 *
 * This file specifies the setup for running end-to-end tests with Cypress.
 *
 * Base URL is set to 'http://localhost:8902' to match the Flask application port.
 * The spec pattern is set to locate test files in the 'tests/cypress/e2e' directory.
 *
 * To run Cypress with this configuration, use the following command:
 * cd tests/cypress
 * npm run cypress:open
 * or
 * npm test
 */


const { defineConfig } = require('cypress')

module.exports = defineConfig({
  video: true,
  e2e: {
    setupNodeEvents(on, config) {},
    baseUrl: 'http://localhost:8902',
    "supportFile": false,
    specPattern: 'e2e/**/*.cy.{js,jsx,ts,tsx}',
  },
})
