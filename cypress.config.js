const { defineConfig } = require('cypress')

module.exports = defineConfig({
  e2e: {
    setupNodeEvents(on, config) {},
    baseUrl: 'http://localhost:8901',
    "supportFile": false,
    specPattern: 'tests/cypress/e2e/**/*.cy.{js,jsx,ts,tsx}',
  },
})
