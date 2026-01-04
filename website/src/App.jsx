import { useState, useEffect } from 'react'
import WitnessList from './components/WitnessList'
import WitnessViewer from './components/WitnessViewer'
import './App.css'

function App() {
  const [witnesses, setWitnesses] = useState([])
  const [selected, setSelected] = useState(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)

  useEffect(() => {
    // Load from static JSONL file built at build time
    fetch('/data/witnesses.jsonl')
      .then(res => {
        if (!res.ok) throw new Error(`Failed to load witnesses: ${res.status}`)
        return res.text()
      })
      .then(text => {
        const polyforms = text
          .trim()
          .split('\n')
          .filter(line => line.length > 0)
          .map(line => JSON.parse(line))
        setWitnesses(polyforms)
        setLoading(false)
      })
      .catch(err => {
        setError(err.message)
        setLoading(false)
      })
  }, [])

  if (loading) {
    return <div className="app loading">Loading witnesses...</div>
  }

  if (error) {
    return <div className="app error">Error: {error}</div>
  }

  return (
    <div className="app">
      <header>
        <div className="header-content">
          <div>
            <h1>Heesch Witness Browser</h1>
            <p className="subtitle">
              Explore polyform tilings and their Heesch numbers
            </p>
          </div>
        </div>
      </header>

      <main>
        <WitnessList
          witnesses={witnesses}
          selected={selected}
          onSelect={setSelected}
        />

        {selected && (
          <WitnessViewer
            witness={selected}
            onClose={() => setSelected(null)}
          />
        )}
      </main>
    </div>
  )
}

export default App
