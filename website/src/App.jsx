import { useState, useEffect } from 'react'
import WitnessList from './components/WitnessList'
import WitnessViewer from './components/WitnessViewer'
import './App.css'

function App() {
  const [witnesses, setWitnesses] = useState([])
  const [selected, setSelected] = useState(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)
  const [refreshing, setRefreshing] = useState(false)

  const loadFromSource = async (url, isJSONL = true) => {
    const res = await fetch(url)
    if (!res.ok) throw new Error(`Failed to load from ${url}`)

    if (isJSONL) {
      // Parse JSONL format (one JSON object per line)
      const text = await res.text()
      return text
        .trim()
        .split('\n')
        .filter(line => line.length > 0)
        .map(line => JSON.parse(line))
    } else {
      // Parse JSON array directly
      return await res.json()
    }
  }

  useEffect(() => {
    // Load from static JSONL file built at build time
    loadFromSource('/data/witnesses.jsonl', true)
      .then(polyforms => {
        setWitnesses(polyforms)
        setLoading(false)
      })
      .catch(err => {
        setError(err.message)
        setLoading(false)
      })
  }, [])

  const handleRefresh = async () => {
    setRefreshing(true)
    setError(null)
    try {
      const polyforms = await loadFromSource(
        'https://hloper--heesch-renderings-web.modal.run/list_full',
        false
      )
      setWitnesses(polyforms)
    } catch (err) {
      setError(err.message)
    } finally {
      setRefreshing(false)
    }
  }

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
          <button
            className="refresh-button"
            onClick={handleRefresh}
            disabled={refreshing}
          >
            {refreshing ? 'Refreshing...' : 'Refresh from Modal'}
          </button>
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
