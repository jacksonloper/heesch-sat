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
    fetch('https://hloper--heesch-renderings-web.modal.run/list_full')
      .then(res => {
        if (!res.ok) throw new Error('Failed to load witnesses')
        return res.json()
      })
      .then(data => {
        setWitnesses(data.polyforms || [])
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
        <h1>Heesch Witness Browser</h1>
        <p className="subtitle">
          Explore polyform tilings and their Heesch numbers
        </p>
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
